% Plot all parameters and cost

figure(2)
for jp1=1:2;
    for jp2=1:9;
        subplot(9,2,jp2+(jp1-1)*9);
            plot(params(jp2+(jp1-1)*9,1:j-1));
            title(par_name(jp2+(jp1-1)*9));
            %hold on
            %plot(ave(jp2+(jp1-1)*9,1:j-1),'r')
    end
end

