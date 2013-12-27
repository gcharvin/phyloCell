function bitPlanes = bitPlanes(im, low, high)
	bitPlanes = bitand(bitshift(im, -low), bitshift(1, (high - low)) - 1);
end
