function text_out = get_asterisks(p,alpha_adjust)
    if p < 0.001/alpha_adjust
        text_out = sprintf('***');
    elseif p < 0.01/alpha_adjust
        text_out = sprintf('**');
    elseif p < 0.05/alpha_adjust
        text_out = sprintf('*');
    else
        text_out = 'ns';
    end
end