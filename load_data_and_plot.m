% reading file and plot R

date = "20-Nov-2023"; % "dd-Mmm-yyyy"
wise = "nonpairwise"; % pairwise / nonpairwise
startat = "sync";  % splay / sync

epsilon = 1e-3;  % overall coupling strength (also K)
epstring = num2str(epsilon);
array = 1;   % array of parameter sets to load

alpharange = "0125"; 

figure 
hold on;

plot_log = 1; % for 1: plot log(1-R), for 0: plot R

clear results

nopar = length(array);
results(nopar).R = [];

if startat == "splay"
  deg = 4; % degree of order parameter to load
else
  deg = 1;
end

for ind = 1:nopar %1:nopar %20
  parsett = num2str(array(ind));
  namet = date +"_Brusselator_"+wise +"_parset" + ...
    parsett+"_"+startat +"_K"+epstring(3:end)+"_alpha" + alpharange +"_N"+num2str(M);
  Mat = readmatrix("data/"+namet + "_R"+deg+".txt");
  results(ind).R = rmoutliers(Mat(2:end,:),"percentiles",[0.5,99.5]);

end

alphasR = Mat(1,2:end); % assuming we load files with same alphas. First argument 0.

black_grey = [0 0 0; 0.4 0.4 0.4];

[~,lenR] = size(results(1).R);
markers = ['x', 'o','+','s','d','^','<','>','p','h'];

rcs_R1 = zeros(size(alphasR));
pol2 = zeros(size(alphasR));
for ind = 1:length(array)
  for ind2 = 2:lenR
    log1_R1 = log(1-results(ind).R(:,ind2));
    Polf = polyfit(results(ind).R(:,1),log1_R1,1);
    if plot_log
      plot(results(ind).R(:,1), log1_R1, "o", 'Color',getcolours(ind2-1)) % 
    else
      plot(results(ind).R(:,1),results(ind).R(:,ind2),'x','linewidth',3,'Color',getcolours(ind2));
    end
    rcs_R1(ind2-1) = Polf(1);
    pol2(ind2-1) = Polf(2);
    hold on;  
  end
end
%%
if plot_log
  for ind = 1:length(array)
    for ind2 = 2:lenR
      plot(results(ind).R(:,1),rcs_R1(ind2-1)*results(ind).R(:,1) + pol2(ind2-1),...
         'Linewidth',3,'Color',getcolours(ind2-1))
      plot(results(ind).R(:,1),rcs_R1(ind2-1)*results(ind).R(:,1) + pol2(ind2-1),":",...
        'Linewidth',3,'Color',[0 0 0])
    end
  end
end

%%
xlabel("\epsilon t");

if plot_log
  ylabel("$\log(1-R_"+num2str(deg)+")$",'Interpreter','latex');
else
  ylabel("$R_"+deg+"$",'Interpreter','latex');
end

if startat == "splay"
  title("Splay bifurcation (set "+parsett+")")
else
  title("In-phase bifurcation (set "+parsett+")")
end

stringalphas = string(num2cell(alphasR));
legendCell = strcat("\Delta\tau = " + stringalphas);
lgd=  legend(legendCell,'NumColumns',4);

%%
figure
hold on;
plot(alphasR,rcs_R1,"LineWidth",3,'Color',getcolours(2))
yline(0)
xlabel("\Delta\tau")
ylabel("slope of $\log(1-R_"+num2str(deg)+")$",'Interpreter','latex');