function filtNs = filterSpikes(Y, dt)
    order=2;
    Ny=1/dt*0.5;
    coF=4;
    coF1=coF/Ny;
    disp(coF1)
    [b c]=butter(order, coF1,'low');
    filtNs=filtfilt(b, c, Y);
