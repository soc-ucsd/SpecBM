function CheckSparsityPattern(Paras)
    [row,col]                = find(Paras.At_sdp);
    Paras.At_sdp_nonzero.row = row;
    Paras.At_sdp_nonzero.col = col;
end