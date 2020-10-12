function [x, y, src] = conncomp2xy( varargin )
%%
if nargin > 0
%     save( strcat( mfilename, '.mat' ) );
else
    load( strcat( mfilename, '.mat' ) );
end
rem_varargin = varargin;
CC = rem_varargin{1};
rem_varargin(1) = [];

tf = cellfun( @islogical, rem_varargin );
ind = 1:CC.NumObjects;
if any( tf )
    ind = find( rem_varargin{tf} );
    rem_varargin(tf) = [];
else
    tf = cellfun( @isnumeric, rem_varargin );
    if any( tf )
        ind = rem_varargin{tf};
        rem_varargin(tf) = [];
    end
end

tf = cellfun( @ischar, rem_varargin );
if any( tf )
    f = find( tf, 1 );
    field = rem_varargin{f};
    rem_varargin(f) = [];
else
    field = 'PixelIdxList';
end


tf = strncmpi( rem_varargin, 'triangle', 3 );
if any(tf)
    dx = 2*[nan; -1; 0; 1; -1; nan];
    dy = 2*[nan; -1; 1; -1; -1; nan]*sqrt(3)/2;    
else %square
    [dy, dx] = ndgrid( [-1 1] );
    dy = [nan; -1; 1; 1; -1; -1; nan];
    dx = [nan; -1; -1; 1; 1; -1; nan];    
end

tf = strcmp( rem_varargin, {'InteriorCentroid'} );
if any(tf)
    rem_varargin(tf) = [];
else
    dx = nan;
    dy = nan;
end

total = sum( cellfun( @numel, CC.(field)(ind) ) );
total = total + numel(dx)*numel(ind);
[y, x, src] = deal( nan(total, 1 ) );
range = 0;
for i=1:numel( ind )
    thing = CC.(field){ind(i)};
    range = range(end) + [1:numel(thing)];
    [y(range), x(range)] = ind2sub( CC.ImageSize, thing );
    src(range) = ind(i);
    %%
    x(range(end)+[1:numel(dx)]) = CC.InteriorCentroid(1, ind(i) ) + dx;
    y(range(end)+[1:numel(dx)]) = CC.InteriorCentroid(2, ind(i) ) + dy;
    %%
    src(range(end)+[1:numel(dx)]) = ind(i);
    range = range(end)+numel(dx);
end

if nargout == 0
    figure( sum( mfilename ) );
    set( clf, 'name', mfilename );
    plot( x, y, '-r', 'ButtonDownFcn', @ButtonDownFcn, 'UserData', src );
end

function ButtonDownFcn( data, event )
%%
[v, m] = min( ( data.XData - event.IntersectionPoint(1) ).^2 + ( data.YData - event.IntersectionPoint(2) ).^2 );
fprintf( '%i\n', data.UserData(m) );

