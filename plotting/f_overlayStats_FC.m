function f_overlayStats_FC(h)

dim = size(h);

hold on;
for hI = 1:dim(1)
    for wI = 1:dim(2)
        if h(hI,wI)
            plot(hI,wI,'w+','MarkerSize',10,'MarkerEdgeColor',[1 1 1],'LineWidth',3);
        end
    end
end

end