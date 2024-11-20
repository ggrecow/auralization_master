%% function  pressure to spl

function result = p2spl(x)
    result = 20.*log10(x./2e-5);
end