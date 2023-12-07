function [Hint, xmin, xmax, ymin, ymax] = load_isochrones(which_oscillator)
    if which_oscillator == "Brusselator"
      Fxx = load( 'BrusselatorFourierAverages.dat');
    elseif which_oscillator == "FitzHughNagumo"
      Fxx = load( 'FHNFourierAverages.dat');
    end

    xxx = Fxx(1,2:end);yxx = Fxx(2:end,1);vxx = Fxx(2:end,2:end)+pi;
    [Xxx,Yxx] = meshgrid(xxx,yxx);Xxx=Xxx'; Yxx=Yxx';vxx=vxx';
    
    xmin = xxx(1); 
    xmax = xxx(end); 
    ymin = yxx(1); 
    ymax = yxx(end);
    
    Hint = griddedInterpolant(Xxx, Yxx, vxx);
end