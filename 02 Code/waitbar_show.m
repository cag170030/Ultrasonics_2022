function waitbar_show(count,L,time,W_bar)
%Input your loop indicator 'count'
%Input your looping seires 'L'
%Input your single loop 'time'

    if count==1
        %W_bar = waitbar(0,'Please wait...');
    waitbar((count/length(L)),W_bar,strcat('Calculating',{' '},num2str(count/length(L)*100),'%'));
    else
        waitbar((count/length(L)),W_bar,strcat('Calculating',{' '},num2str(count/length(L)*100),'%',{' '},'Time remaining',...
            {' '},num2str(floor((length(L)+1-count)*time/60)),'mins',{' '},num2str(((length(L)+1-count)*time/60 ...
            -floor((length(L)+1-count)*time/60))*60),{' '},'Seconds'));
    end
end

