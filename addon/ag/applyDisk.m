function m = applyDisk(m, i, j, r, f)
    [h w] = size(m);
    disk = strel('disk', r, 0).getnhood;
    upward = max(-r, 1 - i);
    downward = min(r, h - i);
    leftward = max(-r, 1 - j);
    rightward = min(r, w - j);
    disk = disk((r + 1 + upward):(r + 1 + downward), (r + 1 + leftward):(r + 1 + rightward));
    m((i + upward):(i + downward), (j + leftward):(j + rightward)) = f(m((i + upward):(i + downward), (j + leftward):(j + rightward)), disk .* r);
end