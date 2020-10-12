function varargout = labelmatrix_sak( varargin )
CC = varargin{1};
varargin = varargin(2:end);
%%
if numel( varargin ) > 0 && ischar( varargin{end} ) && isfield( CC, varargin{end} )
    CC.PixelIdxList = CC.(varargin{end});
    varargin(end) = [];
end

good = [];
tf = cellfun( @islogical, varargin );
if any( tf )
    good = varargin{tf};
    varargin(tf) = [];
end
tf = cellfun( @isnumeric, varargin );
if any( tf )
    good = ismember( 1:CC.NumObjects, varargin{1} );
    varargin(tf) = [];
end

if ~isempty( good )
    CC.PixelIdxList(~good) = {[]};
end

L = labelmatrix( CC );
if nargout > 1
    tmp = vertcat( CC.PixelIdxList{:} );
    [r, c] = ind2sub( CC.ImageSize, tmp );
    axlim = [min(c) max(c) min(r) max(r)];
    L = L(axlim(3):axlim(4), axlim(1):axlim(2) );
    varargout = {L, axlim};
else
    varargout = {L};
end

if nargout == 0
    figure( sum( mfilename ) );
    set( clf, 'Name', mfilename );
    imagesc( L );
    varargout = {};
end