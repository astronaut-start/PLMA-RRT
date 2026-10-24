function [pg, gb] = update_goalbias_pg(pg, gb, success, p)
failNow = ~success;

if gb.count < gb.W
    gb.buf(gb.ptr) = failNow;
    gb.count = gb.count + 1;
    if failNow, gb.fail = gb.fail + 1; end
else
    old = gb.buf(gb.ptr);
    if old, gb.fail = gb.fail - 1; end
    gb.buf(gb.ptr) = failNow;
    if failNow, gb.fail = gb.fail + 1; end
end

gb.ptr = gb.ptr + 1;
if gb.ptr > gb.W, gb.ptr = 1; end

N = min(gb.count, gb.W);
failRate = gb.fail / max(N,1);

if (N >= gb.W) && (failRate > 0.5)
    excess = failRate - 0.5;
    shrink = exp(-p.pg_drop_k * excess);
    pg = max(p.pg_min, pg * shrink);
else
    pg = pg + p.pg_recover_rate * (p.pg_base - pg);
    pg = min(p.pg_max, max(p.pg_min, pg));
end
end
