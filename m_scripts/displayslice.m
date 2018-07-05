function h = displayslice( obj, xslice, yslice, zslice, holdfig, alphacolor)

if nargin <5 | holdfig ==0
    hold off;
else
    hold on
end

sz = size(obj);

ind = find( obj ~= 0);
[i1 i2 i3] = ind2sub(size(obj) , ind);

xmin = min(i2);
xmax = max(i2);
ymin = min(i1);
ymax = max(i1);
zmin = min(i3);
zmax = max(i3);

if xslice & xslice(1) < 0  xslice = round(-xslice * sz(2)); end
if yslice & yslice(1) < 0  yslice = round(-yslice * sz(1)); end
if zslice & zslice(1) < 0  zslice = round(-zslice * sz(3)); end

[Nx, Ny, Nz, Nv] = subvolume(obj, [xmin, xmax, ymin, ymax, zmin, zmax]);

%h=slice( angle(obj23sqsup.object), xslice, yslice, zslice);
h = slice( Nx, Ny, Nz, Nv, xslice, yslice, zslice);
set(h, 'EdgeColor', 'none');

if(alphacolor) alpha(h,'color'); end

hold off;