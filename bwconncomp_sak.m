function varargout = bwconncomp_sak( varargin )
% extended bwconncomp for segmentation
%CC = bwconncomp_sak( bw, cut_off_structure, swath_thickness )
%CC = bwconncomp_sak( bw, cut_off_structure )
%CC = bwconncomp_sak( bw, swath_thickness )
%CC = bwconncomp_sak( CC, cut_off_structure, swath_thickness )
%CC = bwconncomp_sak( CC, cut_off_structure )
%CC = bwconncomp_sak( CC, swath_thickness )
%CC = bwconncomp_sak( CC, ..., '-unique' )

if nargin > 0
    %     save( strcat( mfilename, '.mat'), 'varargin' );
else
    load( strcat( mfilename, '.mat' ) );
end
%%
cut_off = struct( 'Area', 0, 'Perimeter', 0 );

bw_input = islogical( varargin{1} );
if isa( varargin{1}, 'struct' )
    CC = varargin{1};
elseif bw_input
    CC = bwconncomp( varargin{1} );
else
    L = varargin{1};
    ind = unique( L(L>0)' );
    CC = struct( 'ImageSize', size( L ), 'Connectivity', 8, 'NumObjects', max( ind(:) ), 'PixelIdxList', {cell(1, max(ind(:)))} );
    for i = ind
        CC.PixelIdxList{i} = find( L == i );
    end
end
varargin(1) = [];

tf_unique = strncmpi( varargin, '-u', 2 );
if any( tf_unique )
    tf_unique = true;
    varargin(tf_unique) = [];
end
%%
PixelTolerance = 0;
tf = cellfun( @isstruct, varargin );
if any( tf )
    cut_off = varargin{tf};
    varargin(tf) = [];
end

tf = cellfun( @isnumeric, varargin );
if any( tf )
    PixelTolerance = ceil( varargin{tf} / 2);
    varargin(tf) = [];
end

% bw_edge = false( CC.ImageSize );
% tf = strcmp( varargin(1:end-1), 'edge' );
% if any( tf )
%     varargin(tf) = [];
%     bw_edge = varargin{tf};
% end
% bw_edge(1:end, [1 end] ) = true;
% bw_edge([1 end], 1:end ) = true;
% find_edge = find(bw_edge);

[CC.PixelIdxFirst, CC.PixelIdxLast, CC.level] = deal( ones( size( CC.PixelIdxList ), 'uint32' ) );
[CC.Area, CC.Perimeter, CC.Height, CC.MeanRadius, CC.StdRadius, CC.MaxExtent, CC.Neck, CC.AspectRatio, CC.Orientation, CC.Width, CC.Length] = deal( nan( size( CC.PixelIdxList ), 'single' ) );
[CC.Centroid, CC.InteriorCentroid] = deal( zeros( size( CC.PixelIdxList )+[1 0], 'single' ) );
[CC.Radius, CC.PixelIdxListSwath, CC.PixelIdxListContour, CC.rows, CC.cols] = deal( cell( size( CC.PixelIdxList ) ) );

if numel( varargin ) > 0 && isnumeric( varargin{1} ) && isequal( size( varargin{1} ), CC.ImageSize )
    CData = varargin{1};
    varargin(1) = [];
    [CC.MeanIntensity, CC.MaxIntensity, CC.MinIntensity, CC.StdIntensity, CC.MeanEdgeIntensity, CC.MaxEdgeIntensity, CC.MinEdgeIntensity, CC.StdEdgeIntensity] = ...
        deal( nan( size( CC.PixelIdxList ), 'single' ) );
else
    CData = [];
end

% CC.exterior = false( 1, CC.NumObjects );
if ~bw_input
    CC.level = zeros( 1, CC.NumObjects, 'uint16' );
end
if PixelTolerance > 0
    CC.PixelIdxListSwath = cell( size( CC.PixelIdxList ) );
    [rSwath, cSwath] = find( bwmorph( true( 2*PixelTolerance+1 ), 'remove' ) );
    [rSwath, cSwath] = deal( rSwath' - mean(rSwath), cSwath' - mean(cSwath) );
end
%%
level = ones( CC.ImageSize, 'uint32' );
count = 0;

for i=1:CC.NumObjects
    count = count + 1;
    CC.PixelIdxList{count} = uint32( CC.PixelIdxList{i} );
    [CC.rows{count}, CC.cols{count}] = ind2sub( CC.ImageSize, CC.PixelIdxList{i} );
    clim = [min( CC.cols{count} ), max( CC.cols{count} )];
    rlim = [min( CC.rows{count} ), max( CC.rows{count} )];
    %%
    bw = false( [diff(rlim), diff(clim)] + 1 );
    bw( sub2ind( size(bw), CC.rows{count} + 1 - rlim(1), CC.cols{count} + 1 - clim(1) ) ) = true;
    B = builtin("_bwtraceboundarymex", bw, [CC.rows{count}(end) + 1 - rlim(1), CC.cols{count}(end) + 1 - clim(1)], 8 , 0, numel( CC.PixelIdxList{i} ), 0);
    CC.Area(count) = numel( CC.PixelIdxList{i} );
    CC.Perimeter(count) = sum( sqrt( sum( diff(B).^2, 2 ) ) );
    if CC.Perimeter(count) > cut_off.Perimeter && numel( CC.PixelIdxList{i} ) > cut_off.Area
        %     cla; imagesc( bw ); hold on; plot( B(:, 2 ), B(:, 1 ), 'or' );
        CC.PixelIdxListContour{count} = uint32( sub2ind( CC.ImageSize, B(:,1) + rlim(1) - 1, B(:,2) + clim(1) - 1 ) );
        CC.PixelIdxFirst(count) = CC.PixelIdxList{i}(1);
        CC.PixelIdxLast(count) = CC.PixelIdxList{i}(end);
        %%
        hyp2 = ( B(:,1) - B(:,1)').^2 + ( B(:, 2) - B(:, 2)').^2;
        CC.MaxExtent(count) = sqrt( max(hyp2(:)) );
        bw = true( size( bw ) + 2 );
        bw( sub2ind( size( bw ), CC.rows{count} - rlim(1) + 2, CC.cols{count} - clim(1) + 2 ) ) = false;
        d2 = builtin('_bwdistComputeEDT', bw);
        CC.Height(count) = sqrt( max(d2(:)) );
        CC.Centroid(:, count ) = mean( [CC.cols{count}, CC.rows{count}] );
        CC.Radius{count} = sqrt( sum( ( [CC.cols{count}, CC.rows{count}] - CC.Centroid(:, count )' ).^2, 2 ) );
        CC.MeanRadius(count ) = mean( CC.Radius{count} );
        CC.StdRadius(count ) = std( CC.Radius{count} );
        %%
        if size( hyp2, 1 ) > 4
            h_ = hyp2(2:end-1, 2:end-1);
            local_min = h_ < hyp2(2:end-1, 1:end-2 )& h_ < hyp2(2:end-1, 3:end ); % & hyp2(2:end-1, 2:end-1 ) < hyp2(1:end-2, 2:end-1 ) & hyp2(2:end-1, 2:end-1) < hyp2(3:end, 2:end-1 ) ;
            local_min(1:size( local_min, 1 )+1:end) = false;
            local_min(2:size( local_min, 1 )+1:end) = false;
            local_min(1+size( local_min, 1 ) :size( local_min, 1 )+1:end) = false;
            local_min(3:size( local_min, 1 )+1:end) = false;
            local_min(2+size( local_min, 1 ) :size( local_min, 1 )+1:end) = false;
            if any( local_min(:) )
                CC.Neck( count ) = sqrt( min( h_( local_min ) ) );
            end
        end
        %%
        xx = sum( double( CC.cols{i} - CC.Centroid(1, count ) ).^2 );
        yy = sum( double( CC.rows{i} - CC.Centroid(2, count ) ).^2 );
        xy = sum( double( CC.rows{i} - CC.Centroid(2, count ) ).*double( CC.cols{i} - CC.Centroid(1, count ) ) );
        
        [V,D] = eig([ xx xy; xy yy] );
        if D(1) > D(4)
            CC.AspectRatio(count) = sqrt( D(1)/D(4) );
            MajorAxis = V(:,1);
        else
            CC.AspectRatio(count) = sqrt( D(4)/D(1) );
            MajorAxis = V(:,2);
        end
        
        %%
        CC.Orientation( count ) = atan2( MajorAxis(2), MajorAxis(1) );
        tmp = [CC.cols{i} CC.rows{i}]*MajorAxis;
        CC.Length(count) = max(tmp) - min(tmp);
        tmp = [CC.cols{i} CC.rows{i}]*[MajorAxis(2), -MajorAxis(1)]';
        CC.Width(count) = max(tmp) - min(tmp);
        %%
        CC.InteriorCentroid(:, count) = CC.Centroid( :, count );
        if bw( round( CC.Centroid(2, count) ) - rlim(1) + 1, round( CC.Centroid(1, count) ) - clim(1) + 1 )
            %%
            d2 = ( B(:,1) + rlim(1) - 1 - CC.Centroid( 2, count ) ).^2 + ( B(:,2) + clim(1) - 1 - CC.Centroid( 1, count ) ).^2;
            [v2, ind] = min( d2 );
            if v2 > 2
                CC.InteriorCentroid(:, count) = [ B(ind,2) + clim(1),  B(ind,1) + rlim(1)] - 1 ;
            end
        end
        %%
        if PixelTolerance > 0
            [y, x] = deal( B(:,1) - 1 + rlim(1) + rSwath, B(:,2) - 1 + clim(1) + cSwath );
            good = y > 0 & x > 0 & y <= CC.ImageSize(1) & x <= CC.ImageSize(2);
            CC.PixelIdxListSwath{count} = uint32( unique( sub2ind( CC.ImageSize, y(good), x(good) ) ) );
        end
        if ~bw_input
            [CC.level(count), level( CC.PixelIdxList{i} )] = deal( max( level( CC.PixelIdxList{i} ) ) + 1 );
        end
        
        if ~isempty( CData )
            9;
        end
    else
        CC.rows{count} = [];
        CC.cols{count} = [];
        CC.Perimeter(count) = 0;
        if bw_input || tf_unique
            count = count - 1;
        end
    end
    %%
end

if count < CC.NumObjects
    %%
    CC.NumObjects = count;
    fields = setdiff( fieldnames( CC ), {'ImageSize', 'Connectivity', 'NumObjects'} );
    for i = 1:numel(fields)
        CC.(fields{i}) = CC.(fields{i})(:, 1:count );
    end
end

%%
if nargout == 0 || any( strncmpi( varargin, 'plot', 1 ) )
    %%
    fields = fieldnames( CC );
    features = fields( cellfun( @(f) isnumeric( CC.(f) ) && size( CC.(f), 2 ) == CC.NumObjects && size( CC.(f), 1 ) == 1, fields ) );
    %%
    fig = figure( sum( mfilename ) );
    clf;
    set( fig, 'Name', mfilename, 'KeyPressFcn', @KeyPressFcn, 'UserData', struct( 'CC', CC, 'features', {features}, 'CData', CData, 'ind', [] ) );
    set( pan( fig ), 'ActionPostCallback', @AccessKeyPressFcn );
    set( zoom( fig ), 'ActionPostCallback', @AccessKeyPressFcn );
    
    ax.image = findobj( fig, 'Type', 'axes', 'Tag', 'image' );
    if isempty( ax.image )
        ax.image = subplot( 2, 1, 1, 'NextPlot', 'add', 'DataAspectRatio', [1 1 1], 'Tag', 'image', 'ButtonDownFcn', @ButtonDownFcn );
        %         set( zoom( ax.image ), 'ActionPostCallback', @AccessKeyPressFcn );
    end
    ax.histogram = findobj( fig, 'Type', 'axes', 'Tag', 'histogram' );
    if isempty( ax.histogram )
        ax.histogram = subplot( 2, 1, 2, 'NextPlot', 'add', 'Tag', 'histogram', 'yscale', 'log', 'ButtonDownFcn', @ButtonDownFcn );
        ax.histogram.XLabel.String = 'Area';
        ax.histogram.YLabel.String = 'counts';
        set( ax.histogram.XLabel, 'ButtonDownFcn', @ButtonDownFcn );
    end
    
    cla(ax.image);
    if ~isempty( fig.UserData.CData )
        imagesc( ax.image, fig.UserData.CData );
    else
        imagesc( ax.image, labelmatrix( fig.UserData.CC ), 'ButtonDownFcn', @ButtonDownFcn );
        colormap( ax.image, [1 1 1; .7*ones( fig.UserData.CC.NumObjects, 3 )] );
    end
	axis(ax.image, 'tight' );
    [x, y, ind] = conncomp2xy( fig.UserData.CC, 'PixelIdxListContour' );
    
    %%
    plot( ax.image, x, y, '-k', 'UserData', ind, 'ButtonDownFcn', @ButtonDownFcn, 'DisplayName', 'Contour' );
    plot( ax.image, fig.UserData.CC.Centroid(1,:), fig.UserData.CC.Centroid(2,:), 'mo', 'UserData', 1:fig.UserData.CC.NumObjects, 'ButtonDownFcn', @ButtonDownFcn, 'DisplayName', 'Centroid' );
    tf = any( fig.UserData.CC.Centroid ~= fig.UserData.CC.InteriorCentroid, 1 );
    plot( ax.image, fig.UserData.CC.Centroid(1,tf), fig.UserData.CC.Centroid(2,tf), '*g', 'UserData', 1:fig.UserData.CC.NumObjects, 'ButtonDownFcn', @ButtonDownFcn, 'DisplayName', 'Interior Centroid' );
    quiver( fig.UserData.CC.Centroid(1,:), fig.UserData.CC.Centroid(2,:), cos( fig.UserData.CC.Orientation ), sin( fig.UserData.CC.Orientation ), .2, '-r', 'ButtonDownFcn', @ButtonDownFcn, ...
        'DisplayName', 'Direction', 'Parent', ax.image, 'UserData', 1:fig.UserData.CC.NumObjects, 'Visible', 'off' );
    leg = findobj( ax.image.Parent, 'Type', 'legend', 'UserData', ax.image );
    if isempty( leg )
        leg = legend(ax.image);
        set( leg, 'ItemHitFcn', @LegendItemHitFcn );
    end
    
    update_axes( fig );
    %     update_fig( figure(fig) );
    varargout = { CC, fig };
else
    varargout = { CC };
end

varargout = varargout( 1:nargout );

function KeyPressFcn( fig, event )
switch event.Character
    case '?'
        fprintf( 'z: toggle zoom\n' );
        fprintf( 'Z: toggle pan\n' );
        fprintf( 'f: change feature \n' );
    case  {'f', 'F'}
        %%
        for i = 1:numel(fig.UserData.features)
            fprintf( '%i: %s\n', i, fig.UserData.features{i} );
        end
        i = input( 'enter # for feature: ' );
        ax = findobj( fig, 'Tag', 'histogram' );
        ax.XLabel.String = fig.UserData.features{i};
        update_axes( fig );
    case 'z'
        zoom( gca );
    case 'Z'
        pan( gca );
end
%%

function ButtonDownFcn( data, event )
%%
fig = gcf;
switch data.Type
    case 'image'
        fig.UserData.ind = data.CData( round( event.IntersectionPoint(2) ), round(event.IntersectionPoint(1) ) );
    case 'text'
        tf = strcmp( data.String, fig.UserData.features );
        if event.Button == 1
            tf = tf( [end 1:end-1] );
        else
            tf = tf( [2:end 1] );
        end
        data.String = fig.UserData.features{tf};
    case {'axes', 'histogram'}
        %%
        ax.histogram = findobj( fig, 'type', 'axes', 'tag', 'histogram' );
        ax.image = findobj( fig, 'type', 'axes', 'tag', 'image' );
        hist = findobj( ax.histogram, 'type', 'histogram' );
        
        BinCenters = mean( hist.BinEdges([1:end-1; 2:end] ) );
        [v, m] = min( ( BinCenters - event.IntersectionPoint(1) ).^2 );
        f = find( fig.UserData.CC.( ax.histogram.XLabel.String ) > hist.BinEdges(m) & fig.UserData.CC.( ax.histogram.XLabel.String ) <= hist.BinEdges(m+1) & ...
            fig.UserData.CC.Centroid(1,:) > min( xlim( ax.image ) ) &  fig.UserData.CC.Centroid(1,:) < max( xlim( ax.image ) ) & ...
            fig.UserData.CC.Centroid(2,:) > min( ylim( ax.image ) ) &  fig.UserData.CC.Centroid(2,:) < max( ylim( ax.image ) ) );
        if ~isempty( f )
            [v, i] = min( sum( ( [interp1( hist.BinEdges(m+[0 1]), xlim( ax.image ), event.IntersectionPoint(1) ); ...
                interp1( ylim( ax.histogram ), ylim( ax.image ), event.IntersectionPoint(2) )] - fig.UserData.CC.Centroid(:, f) ).^2, 1 ) );
            fig.UserData.ind = f(i);
        end
    case {'line', 'quiver'}
        [v, m] = min( sum( ( [data.XData; data.YData]' - event.IntersectionPoint([1 2] ) ).^2, 2 ) );
        fig.UserData.ind = data.UserData(m);
end
update_axes( fig )


function update_axes( fig )
ax.histogram = findobj( fig, 'type', 'axes', 'tag', 'histogram' );
hist = findobj( ax.histogram, 'type', 'histogram' );
if isempty( hist ) || ~strcmp( hist.UserData, ax.histogram.XLabel.String )
    cla(ax.histogram);
    hist = histogram( ax.histogram, fig.UserData.CC.( ax.histogram.XLabel.String ), 'ButtonDownFcn', @ButtonDownFcn, 'UserData', ax.histogram.XLabel.String, 'Tag', 'histogram1' );
    axis( ax.histogram, 'tight' );
    ylim( ax.histogram, [min( ylim( ax.histogram ) )*.9, max( ylim( ax.histogram ) )] );
end

if fig.UserData.ind > 0
    fprintf( 'ind %i\n', fig.UserData.ind );
    for i=1:numel(fig.UserData.features)
        fprintf( '\t%s %s\n', fig.UserData.features{i}, num2str( fig.UserData.CC.(fig.UserData.features{i})(fig.UserData.ind) ) );
    end
    %%
    ax.image = findobj( fig, 'type', 'axes', 'tag', 'image' );
    focus = findobj( ax.image, 'Tag', 'focus' );
    if isempty( focus )
        focus = plot( ax.image, nan, nan, 'hk', 'MarkerFaceColor', 'c', 'Tag', 'focus' );
    end
    set( focus, 'XData', fig.UserData.CC.Centroid(1, fig.UserData.ind ), 'YData', fig.UserData.CC.Centroid(2, fig.UserData.ind ), ...
        'DisplayName', sprintf( 'ind %i', fig.UserData.ind ) );
    
    focus = findobj( ax.histogram, 'Tag', 'focus' );
    if isempty( focus )
        focus = plot( ax.histogram, nan, nan, 'hk', 'MarkerFaceColor', 'y', 'Tag', 'focus' );
    end
    BinCenters = mean( hist.BinEdges([1:end-1; 2:end] ) );
    [v, m] = min( ( BinCenters - double( fig.UserData.CC.( ax.histogram.XLabel.String )(fig.UserData.ind) ) ).^2 );
    set( focus, 'XData', interp1( xlim( ax.image ), hist.BinEdges(m+[0 1] ), fig.UserData.CC.Centroid(1,fig.UserData.ind) ), ...
        'YData', interp1( ylim( ax.image ), ylim( ax.histogram ), fig.UserData.CC.Centroid(2,fig.UserData.ind) ) );
end
%%
function LegendItemHitFcn( plt, event )
if strcmp( event.Peer.Visible, 'off' )
    set( event.Peer, 'Visible', 'on' );
else
    set( event.Peer, 'Visible', 'off' )
end

function AccessKeyPressFcn( fig, event)
hManager = uigetmodemanager(fig);
[hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
%%
hist1 = findobj( fig, 'Tag', 'histogram1' );
%%
axs = axis( findobj( fig, 'Tag', 'image' ) );
sub = fig.UserData.CC.Centroid(1,:) > axs(1) & fig.UserData.CC.Centroid(1,:) < axs(2) & fig.UserData.CC.Centroid(2,:) > axs(3) & ...
    fig.UserData.CC.Centroid(2,:) < axs(4);
counts = histcounts( fig.UserData.CC.(hist1.Parent.XLabel.String)(sub), 'BinEdges', hist1.BinEdges );
zoom_histogram = findobj( hist1.Parent, 'Tag', 'zoom_histogram' );
if isempty( zoom_histogram )
    zoom_histogram = plot( hist1.Parent, nan, nan, '+r', 'Tag', 'zoom_histogram' );
end
set( zoom_histogram, 'XData', mean(  hist1.BinEdges( [1:end-1; 2:end] ) ), 'YData', counts );


%%
set(fig, 'KeyPressFcn', @KeyPressFcn );

