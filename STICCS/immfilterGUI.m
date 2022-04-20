function movie = immfilterGUI(movie,filterType,msg,ax)


cla(ax)
ylim(ax,[0,1])
xlim(ax,[0,1])
ph = patch(ax,[0 0 0 0],[0 0 1 1],[0.67578 1 0.18359]); %greenyellow
th = text(ax,1,1,[msg '...0%'],'VerticalAlignment','bottom','HorizontalAlignment','right');


if strcmp(filterType,'Mobile')
    movie = fft(double(movie),[],3);
    avPixInt = mean(mean(movie(:,:,1)))/size(movie,3);
    movie(:,:,2:end) = 0;
    movie = uint16(real(ifft(movie,[],3))+avPixInt);
% 'F' option FFT filters the whole movie at once
% fastest for regular sized movies
elseif strcmp(filterType,'F')
    movie = fft(double(movie),[],3);
    ph.XData = [0 1/3  1/3 0];
    th.String = sprintf([msg '...%.0f%%'],round((1/3)*100));
    drawnow %update graphics
    avPixInt = mean(mean(movie(:,:,1)))/size(movie,3);
    ph.XData = [0 2/3  2/3 0];
    th.String = sprintf([msg '...%.0f%%'],round((2/3)*100));
    drawnow %update graphics
    movie(:,:,1) = 0;
    movie = (real(ifft(movie,[],3))+avPixInt);
    ph.XData = [0 1 1 0];
    th.String = sprintf([msg '...%.0f%%'],round(100));
    drawnow %update graphics
    % 'FM' option FFT filters each pixel location
    % one at a time.  Use for large movies... uses less
    % memory than 'F' option, but probably a little slower
elseif strcmp(filterType,'FM')
    AvIntSum = 0;
    for i = 1:size(movie,1)
        
        for j = 1:size(movie,2)
            fftPixelTrace = fft(double(movie(i,j,:)));
            AvIntSum = AvIntSum + fftPixelTrace(1);
            fftPixelTrace(1) = 0;
            movie(i,j,:) = int16(ifft(fftPixelTrace));
        end
        ph.XData = [0 i/size(movie,1)  i/size(movie,1) 0];
        th.String = sprintf([msg '...%.0f%%'],round(i/size(movie,1)*100));
        drawnow %update graphics
    end
    movie = movie + AvIntSum/numel(movie);
    % If filtertype is a number, than use this as a moving average window size
elseif strcmp(filterType,'F0')
    for i = 1:size(movie,1)
        for j = 1:size(movie,2)
            fftPixelTrace = fft(double(movie(i,j,:)));
            fftPixelTrace(1) = 0;
            movie(i,j,:) = int16(ifft(fftPixelTrace));
        end
        ph.XData = [0 i/size(movie,1)  i/size(movie,1) 0];
        th.String = sprintf([msg '...%.0f%%'],round(i/size(movie,1)*100));
        drawnow %update graphics
    end
    % If filtertype is a number, than use this as a moving average window
    % size
elseif isnumeric(filterType)
    % This old way is commented out
    % old way #1 -- using matlab's filter
    %     movieMean = filter(ones(1,filterType)/filterType,1,double(movie),[],3);
    %     movie = int16(double(movie) - movieMean);

    % old way number #2 -- using moving_average
    F = round((filterType - 1)/2);
    for i = 1:size(movie,1)
        for j = 1:size(movie,2)
            averagePixelTrace = (moving_average(double(movie(i,j,:)),F));
            movie(i,j,:) = movie(i,j,:)-averagePixelTrace;
        end
        ph.XData = [0 i/size(movie,1)  i/size(movie,1) 0];
        th.String = sprintf([msg '...%.0f%%'],round(i/size(movie,1)*100));
        drawnow %update graphics
    end
elseif strcmp(filterType,'none')
    % do nothing
else
    error('Filtertype must be ''F'', ''FM'', ''none'', or for a moving average, a number.');
end


