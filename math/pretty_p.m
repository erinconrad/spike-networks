function pp = pretty_p(p,adjust_alpha)
    ast = get_asterisks(p,adjust_alpha);
    
    if p < 0.001
        pp = sprintf('p < 0.001%s',ast);
    elseif p < 0.01
        pp = sprintf('p = %1.3f%s',p,ast);
    elseif p < 0.05
        pp = sprintf('p = %1.3f%s',p,ast);
    else
        pp = sprintf('p = %1.2f%s',p,ast);
    end
end