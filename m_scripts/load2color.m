function im = load2color(filename)

im = double(h5read(filename, '/entry/data/data/'));
im = im(:,:,1)-im(:,:,2);


