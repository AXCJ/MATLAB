function plothg(inputX, inputY)
h = plot(0,0,'XDataSource','x','YDataSource','y');

for k=1:length(inputX)
    x=inputX(1:k);
    y=inputY(1:k);
    refreshdata(h)
    pause(1e-9)
    drawnow
end

end

