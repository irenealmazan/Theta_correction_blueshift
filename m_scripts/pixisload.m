function ccd = pixisload(number, displayflag, range)

if nargin<2
    displayflag=1;
end

if nargin<3
    range=[1 1024 1 1024];
end

ccd = hdfread(['image' num2str(number, '%1.5d') '.hdf'], '/entry1/data/data', 'Index', {[1 1],[1 1],[1024 1024]});
ccd = ccd(range(3):range(4), range(1):range(2));

if(displayflag)
    imagesc(ccd)
    axis equal tight
end