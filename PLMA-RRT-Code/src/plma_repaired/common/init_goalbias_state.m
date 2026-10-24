function gb = init_goalbias_state(W)
gb.W     = W;
gb.buf   = false(W,1);
gb.ptr   = 1;
gb.count = 0;
gb.fail  = 0;
end
