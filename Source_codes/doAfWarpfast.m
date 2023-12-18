
function [pts, tree, im, af] = doAfWarpfast(pts, tree, im, fiximsz)

imsz = size(im);

imold = im;
im = zeros(max(fiximsz(1), imsz(1)), max(fiximsz(2), imsz(2)));
im(1:imsz(1), 1:imsz(2)) = imold;
clear imold;

posf = pts.fix;
posm = pts.mov;

posm = [posm(:, [2 1]) ones(size(posm, 1), 1)];
posf = [posf(:, [2 1]) ones(size(posf, 1), 1)];

mid = [size(im)./2 1];

posm(:,1) = posm(:,1) - mid(1);
posm(:,2) = posm(:,2) - mid(2);
% posm(:,3) = posm(:,3) - mid(3);

af = posm\posf;

af(3,1) = af(3,1) - mid(1);
af(3,2) = af(3,2) - mid(2);
% af(4,3) = af(4,3) - mid(3);

im = affine_transform(im, inv(af'),1);

padx = round(fiximsz(2) - imsz(2));
pady = round(fiximsz(1) - imsz(1));
% padz = round(fiximsz(3) - imsz(3));

dxn = round(padx.*(padx<0));
dyn = round(pady.*(pady<0));
% dzn = round(padz.*(padz<0));

im = im(1:end+dyn, 1:end+dxn);

posf = pts.fix;
posm = pts.mov;

posm = [posm ones(size(posm, 1), 1)];
posf = [posf ones(size(posf, 1), 1)];

af = posm\posf;

af(abs(af)<1e-9) = 0;
af(abs(af-1)<1e-9) = 1;

a = [pts.mov, ones(size(pts.cur,1), 1)];
a = a*af;
pts.cur = a(:,1:2);

if ~isempty(tree.mov)
    a = [tree.mov, ones(size(tree.mov,1), 1)];
    a = a*af;
    tree.cur = a(:,1:2);
end

