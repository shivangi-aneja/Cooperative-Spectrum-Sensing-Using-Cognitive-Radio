function y=MyProj2()
    clear all
    close all
    clc
    %----------- frequency range to be scanned for a signal ------------%
    range_start=54* 10^6;        %frequency in Hz i.e. 54 MHz
    range_end= 698 * 10^6;     %end range of 698 MHz
	
    %------------Primary User information-------------%
    PX=rand*50+30;    %primary user's location
    PY=rand*100;
    
    Fs= 1400*10^6;    %sampling frequency 1400MHz
    
    Primary_User_Freq= (range_start + (range_end-range_start)*rand(1));         % Carrier frequency- some random singnal in the range which has to be detected

    t = 0:0.00001:0.01;           
    msg= 10*randn(1,size(t,2));                                                 % real valued gaussian sginal
    Primary_User_Signal= ammod(msg,Primary_User_Freq,Fs);                       %final primary user signal
	
    %-------------------- periodogram of original signal ------------------- %
    figure(1)
    Pxx=periodogram(Primary_User_Signal);
    Hpsd = dspdata.psd(Pxx,'Fs',Fs);
    plot(Hpsd);
    grid on;
    title('Original Signal - Test');    
   
    %------------Secondary users(nodes) --------------%
    NoOfNodes=10;                                                                 % number of secondary users
    nodes = 1:NoOfNodes;

    fusionNode=floor(((rand+1)*NoOfNodes)/2);                                     % the fusion node at which all data is collected at the end

    X=rand(size(nodes)).*100;                                                     % random locations of nodes
    Y=rand(size(nodes)).*100;
    showAllNodes(NoOfNodes,X,Y,PX,PY,fusionNode)

    velocityi= floor(rand(size(nodes)).*3);                                       % random velocities of the nodes
    velocityj= floor(rand(size(nodes)).*3);
    
    Pf = mod(rand(size(nodes)),0.1);      %Velocity of the fusion node

    %---------- getting data from all the secondary users and displaying ---------------%
    disp('Running Simulation');
    x=0;
    run=0;
    allfinal=zeros(1,86);
    for time=1:500
    
	%doubttt
    if(mod(time,50)==0)
        run=run+1;
        [new_decision,new_energy]=secondary_users(NoOfNodes,X,Y,PX,PY,Primary_User_Signal,range_start,range_end,Fs,velocityi,velocityj,Pf,run); %#ok<*NASGU>
        
        % [new_decision,new_energy] - this is the local decision of each of the secondary users.
        
        %fusion centre code - here it starts!
        %computing weignts using energies
        p=0;
        k=0;
        sum=zeros(86);
        %final_decision=zeros(86);
        %final_decision=zeros(1,86);
        for k=1:86      %for each channel
            
           for p=2:(NoOfNodes +1)
            sum(k)=sum(k)+ new_energy(p,k);    %take sum of the energies of all the secondary users in a particular channel
           end
           
           avrg(k)=sum(k)/NoOfNodes;    %average of the energies of all the secondary users in a particular channel
           
           temp =0;
           
           for j=2:(NoOfNodes+1)
              weight(j)= new_energy(j,k)/avrg(k);      %weight of a particular secondary user
              final_weight(j,k)=weight(j);     %final weight of a secondary user in the kth channel
              temp = temp + weight(j)* new_decision(j,k);   %The weight assigned to every secondary user is multiplied to the local decision value and the cumulative sum obtained from all the secondary users is used to determine the final decision of the fusion centre.
           end
           
           if(temp>0)     %computing the final decision based on the cumulative  value of temp
               final_decision(k) = 1;
           else
               final_decision(k) = 0;
           end
           
        end      %end of the channels loop (k-> 1:86)
       
        allfinal = vertcat(allfinal,final_decision);
        save('Final','allfinal');
        %finalmartix=vertcat(finalmatrix, final_decision);
        %disp(final_decision);
        %disp(final_decision);
        %fusion centre code ends
        %disp(decision);
        %disp(energy);
        % --- The above function call generates periodograms for all the
        %         secondary users and takes decisions based on calculated threshold
        % --- All the variables regarding a secondary user are stored in the 10
        %       ".mat" files generated. They can be loaded using the "load()" function of MATLAB
        % --- The decision is an array of (698-54)/8 =80 (i took 100)
        %         channels. it consists of zeroes and ones        How??
        % --- The threshold consists of thresholds used to calculate and
        %         distance is the distance of the SU from the primary user
        
        if(mod(time,100)==0)
             x=x+1;
            [velocityi velocityj]=changeVelocities(velocityi,velocityj,x);
        end
        
        % *************** THE code/function call, must be here *********************** %
        [X Y]=changeDistances(velocityi,velocityj,X,Y);
        showAllNodes(NoOfNodes,X,Y,PX,PY,fusionNode);
        pause(1);
    %save('final_matrix','finalmatrix');
    end
    
    end     %this is the end of the main (time, 1:500) loop for all the 10 runs.
    xlswrite('Data.xlsx',allfinal,'Fusion','B2');
    disp('This is reached. This is the end of all the runs in the simulation.');
end

%function for local decisions of each of the secondary users. 
function [decsn,enrgy]=secondary_users(NoOfNodes,X,Y,PX,PY,Primary_User_Signal,range_start,range_end,Fs,velocityi,velocityj,Pf,run)

    
    %------- calculating distance,signals,periodograms at each node ------%
     alldecisions=zeros(1,100);
     allenergies=zeros(1,100);
        
     
    for i=1:NoOfNodes
      
       distance=sqrt( (X(i)-PX)^2 + (Y(i)-PY)^2) ;                       %distance of each node from primary user is placed in the distance array
       snr_at_node=250/(distance);
       signal_at_node=awgn(Primary_User_Signal,snr_at_node);              %adding AWGN (noise) to the signal  y(n)=s(n)+w(n)
       
       %***************** calculating decision metric and threshold *****************%
 
        L=size(Primary_User_Signal,2);
        Threshold= (qfuncinv(Pf(i))./sqrt(L))+ 1;     
       
       %************ Finding periodogram at the secondary user *********%
        periodogram_at_node=periodogram(signal_at_node);                      %calculating the periodograms of all the secondary users
        figure(3);
        Hpsd = dspdata.psd(periodogram_at_node,'Fs',Fs);
        subplot(5,2,i);
        plot(Hpsd);
        grid on;
        title(sprintf('Periodogram at node- %d; dis- %d vel-%di+%dj',i,distance,velocityi(i),velocityj(i)));
        
       %******************* Coming to a decision ******************%
       
        %-----Every secondary user scans the range to find the primary signal----%
        occupied=scanRange(range_start,range_end,periodogram_at_node,Threshold);                         
        Threshold=Threshold*10;
        g=1;
        h=1;
        k=1;
        decision(100)=0;
        power(100)=0;
        count=0;
        channel_width =7;   %MHz
        max=Fs/2;
        min = max/length(occupied);
        width = channel_width/min;
        width=width*10^6;
        width = ceil(width);
%         disp(channel_width);
%         disp(max);
%         disp(min);
%         disp(width);
        while(g<length(occupied))
           sum=0;
           energy=0;
           for h=1:width
               if(g+h<length(occupied))
                    sum=sum+occupied(g+h);
                    energy = energy+periodogram_at_node(g+h);
               end
           end
           power(k)=energy;
           if(sum>=width/2)
               decision(k)=1;
               count=count+1;
           else
               decision(k)=-1;
           end
%            allenergies(i,k)=energy;
%            alldecisions(i,k)= decision(k);
           
           k=k+1;
           g=g+width;
        
        end
         decsn=decision;
         enrgy=power;
         alldecisions =vertcat(alldecisions,decision);
         allenergies=vertcat(allenergies,power);
        %save(sprintf('variables_at_node- %d in run - %d ',i,run),'distance','decision','count','power');
        %save(sprintf('variables_at_node- %d ',i),'distance','decision','count','power');
    end
	
       
     decsn=alldecisions;
     enrgy=allenergies;
    % save(sprintf('this is final at run - %d',run),'alldecisions','allenergies');
    %WRITING TO EXCEL
    str1='B';
    str2='A';
    num=3+ (run-1)*13;
    cell=strcat(str1, int2str(num));
    cell2=strcat(str2,int2str(num-1));
    %disp(cell);
     %disp(cell2);
   runvalue=strcat('Run',int2str(run));
 %   permanent=[runvalue;'Garbage Value';'Node1';'Node2';'Node3';'Node4';'Node5';'Node6';'Node7';'Node8';'Node9';'Node10'];
    xlswrite('Data.xlsx',{runvalue},'Energy',cell2);
    xlswrite('Data.xlsx',allenergies,'Energy',cell);
    xlswrite('Data.xlsx',{runvalue},'LocalDecision',cell2);
    xlswrite('Data.xlsx',alldecisions,'LocalDecision',cell);     %these are for the local decisions of each run!
    extra=strcat('Ended ',runvalue);
    disp(extra);     %here is it EndedRun(n)
    %END WRITING TO EXCEL



end

%------------Plot of the nodes ------------%
function showAllNodes(NoOfNodes,X,Y,PX,PY,fusionNode)
	
    figure(2);
    subplot(1,1,1,'replace');
    axis([-10 120 -10 120]);
    title('All Nodes');
    %plot(X,Y,'.',PX,PY,'*',X(fusionNode),Y(fusionNode),'.');
    grid on;
    hold on;
    plot(PX,PY,'r*');
    plot(X(fusionNode),Y(fusionNode),'o');
    for i=1:NoOfNodes
        plot(X(i),Y(i),'.');
    end
end

function [velocityi velocityj]=changeVelocities(velocityi,velocityj,x)
if(mod(x,4)==0)
    velocityi = (-1*velocityi);
    if (mod(x,2)==0)
        velocityj = (-1*velocityj );
    end
elseif (mod(x,4)==1)
    velocityj = (-1*velocityj );
    if (mod(x,2)==1)
        velocityi = (-1*velocityi );
    end
end

end

function [di dj]=changeDistances(velocityi,velocityj,X,Y)
    di=velocityi;
    dj=velocityj;
    X=X+di;
    Y=Y+dj;
    for j=1:size(X)
        X(j)=X(j)+di(j);
        if(X(j)>=95)
            X(j)=X(j)-2*di(j);
        elseif(X(j)<=5)
            X(j)=X(j)+2*di(j);
        end
        Y(j)=Y(j)+dj(j);
        if(Y(j)>=95)
            Y(j)=Y(j)-2*dj(j);
        elseif(Y(j)<=5)
            Y(j)=Y(j)+2*dj(j);
        end
    end
    di=X;
    dj=Y;
end

function y=scanRange(rangestart,rangeend,Pxx,threshold)

     occupied(length(Pxx))=0;
    i=1; threshold=threshold*10;
%     disp('length of Pxx');
%     disp(length(Pxx));
    while (i<length(Pxx))
        
        if(Pxx(i)>threshold)
            occupied(i)=1;
        end
        i=i+1;
    end

    y=occupied;
    
end



