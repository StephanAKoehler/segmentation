function varargout = L_fragment( varargin )
% CC_frag = L_fragment( L_input )
% CC_frag = L_fragment( CC_input )


if nargin == 0
    %
    bw.input = [
        5     5     0     5     5     5     0     0     0     0     0     0     0     0     0
        5     5     5     5     5     5     0     0     0     0     0     0     0     0     0
        5     5     0     5     5     5     0     0     0     0     0     0     0     0     0
        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
        0     0     0     0     0     0     2     2     0     2     2     0     0     0     0
        0     0     0     0     0     0     2     2     2     2     2     0     0     0     0
        0     0     0     0     0     0     2     2     2     2     2     0     0     0     0
        0     0     0     0     0     0     2     2     0     2     2     0     0     0     0
        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
        0     0     0     0     0     0     0     0     0     0     0     3     0     0     0
        0     0     0     0     0     0     0     0     0     0     0     0     0     4     0
        0     0     0     0     0     0     0     0     0     0     0     0     0     4     0
        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0];
    bw.input = logical( [
        0     0     0     0     0     1     1     0     1     1     1     1     1     1     0     1     1     1     1     1
        0     0     0     0     0     1     1     1     1     1     1     1     1     1     0     1     1     1     1     1
        1     1     0     0     0     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
        1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
        1     1     1     1     1     1     1     1     1     1     1     1     0     1     1     1     1     1     1     1
        1     1     0     0     0     1     1     0     1     1     1     1     0     1     1     1     1     1     1     1
        1     0     0     0     0     1     1     0     1     1     1     1     0     1     1     1     1     1     1     1
        ] );
    
    bw.input = logical( ...
        [0     0     0     0
        0     1     0     0
        0     1     0     0
        0     0     0     0] );
    %         bw_embed = false( size(bw.input) + 2 );
    %         bw_embed(2:end-1, 2:end-1) = bw.input;
    %         bw.input = repmat( bw_embed, 100, 60 );
    %     bw.input = bwlabeln( bw.input );
    varargin = { bw.input, bw.input };
    
    load( strcat( mfilename, '.mat' ) );
    varargin{2} = struct( 'iter', 10000 );
else
%     save( strcat( mfilename, '.mat' ) );
end
cutoff.iter = 1000;
cutoff.Area = 30;
if numel( varargin ) > 1 && isstruct( varargin{2} )
    for fields = intersect( fieldnames( varargin{2} ), fieldnames( cutoff ) )
        cutoff.(char(fields)) = varargin{2}.(char(fields));
    end
end

%%
CC.input = varargin{1};
L.input = labelmatrix( CC.input );
%%
bw.input = (L.input ~= 0 );
bw_ = true( CC.input.ImageSize + 2);
bw_(2:end-1,2:end-1) = ~bw.input;
% D = sqrt( builtin('_bwdistComputeEDT', bw.input) );
d_bwdist = sqrt( builtin('_bwdistComputeEDT', bw_) );
d_bwdist = d_bwdist( 2:end-1, 2:end-1 );
bw.watershed = ( watershed( -imhmin( d_bwdist, 1 ) ) ==  0 ) & bw.input;
%%
% figure(11); clf; imagesc( d_bwdist.*single( bw.input ) ); axis( [ 470 515 285 330] ); colorbar
% %%
% figure(11); clf; imagesc( d_bwdist ); axis( [ 470 515 285 330] ); colorbar
% %%
% figure(11); clf; imagesc( bw.watershed );
%%
bw.pieces = ismember( L.input, L.input( bw.watershed ) ) & ~bw.watershed;
level = uint16( bw.pieces );
CC.pieces = bwconncomp( bw.pieces );
CC.pieces.level = ones( 1, CC.pieces.NumObjects, class(level) );
L.pieces = labelmatrix( CC.pieces );
%%
d = L.pieces( [2:end, end], : );
u = L.pieces( [1 1:end-1], : );
l = L.pieces( :, [2:end, end] );
r = L.pieces( :, [1 1:end-1] );

comp_reg = @( a, b ) a~=b&a~=0&b~=0;

bw.lr = comp_reg( l, r );
bw.ud = comp_reg( u, d );

% imagesc( L.pieces ); hold on; [y, x] = find(  bw.lr | bw.ud ); plot( x, y, 'or' )
pieces12 = [l( bw.lr & bw.watershed), r( bw.lr & bw.watershed ); u( bw.ud & bw.watershed ), d( bw.ud & bw.watershed )];
pieces12_max = max( pieces12(:) );
pieces12_ind = sub2ind( pieces12_max([1 1]), pieces12(:,1), pieces12(:,2) ) + sub2ind( pieces12_max([1 1]), pieces12(:,2), pieces12(:,1) );
L.watershed = zeros( CC.input.ImageSize );
sm = sum( sum(  bw.lr & bw.watershed ) );
L.watershed( bw.lr & bw.watershed ) = pieces12_ind(1:sm );
L.watershed( bw.ud & bw.watershed ) = pieces12_ind(1+sm:end );
% clf;
% imagesc( L.watershed + bw.pieces + bw.input );
% set( gca, 'clim', [0 3] );
%%

cut_length_pieces12 = histcounts2( [pieces12(:,1); pieces12(:,2)], [pieces12(:,2); pieces12(:,1)], ...
    1:CC.pieces.NumObjects + 1, 1:CC.pieces.NumObjects + 1 );
%%
CC.pieces.CutLength = sum( cut_length_pieces12 );
CC.pieces.MaxCutLength = max( cut_length_pieces12 );
%%
pieces12 = [l( bw.lr ), r( bw.lr ); u( bw.ud ), d( bw.ud )];
pieces12 = sortrows( sort( pieces12, 2 ), [1 2] );
tf = [true; any( pieces12(1:end-1, : ) ~= pieces12( 2:end, : ), 2 )];
pieces12 = pieces12( tf, : );
input_pieces = unique( [L.input(L.pieces>0) L.pieces( L.pieces > 0 )], 'rows' );

%%
Neighbors = cell( 1, CC.pieces.NumObjects );
for i = 1:CC.pieces.NumObjects
    Neighbors{i} = pieces12( any( ismember( pieces12, i ), 2 ), : );
end
%%

for i = unique( input_pieces(:,1) )'
    %%
    pieces = input_pieces( input_pieces( :, 1) == i, 2 );
    bnry = get_binary( numel(pieces) );
    bnry = bnry( sum( bnry, 2 ) > 1, : );
    if ~isempty( bnry )
        [r, c] = ind2sub( CC.input.ImageSize, CC.input.PixelIdxList{i} );
        r_lim = [min(r(:)) max( r(:) ) ];
        c_lim = [min(c(:)) max(c(:)) ];
        bw_input = L.input( r_lim(1):r_lim(2), c_lim(1):c_lim(2) ) == i;
        L_pieces = L.pieces( r_lim(1):r_lim(2), c_lim(1):c_lim(2) ).*cast( bw_input, 'like', L.pieces );
        for j = 1:numel( pieces )
            %%
            tmp = pieces';
            tmp(j) = 0;
            tf = bnry(:,j);
            bad = false( size( tf ) );
            bad(tf) = ~any( bnry(tf,:) & ismember( tmp, Neighbors{pieces(j)} ), 2 );
            bnry(bad, : ) = 0;
        end
        %%
        bnry = bnry(any(bnry,2), : );
        if size( bnry, 1 ) < cutoff.iter
            for j = 1:size( bnry, 1 )
                %%
                p = pieces( bnry(j, : ) );
                tf_connected = false( size( p ) );
                tf_connected(1) = true;
                for k = 1:numel( p )
                    lngth = sum( tf_connected );
                    tf_connected = tf_connected | ismember( p, vertcat( Neighbors{p(tf_connected)} ) );
                    if sum( tf_connected ) == lngth
                        break
                    end
                end
                if all( tf_connected ) && k == numel(p)
                    %%
                    tf = ismember( L_pieces, p );
                    %%
                    L_conn = zeros( size(L_pieces)+2, class( L_pieces ) );
                    L_conn(2:end-1, 2:end-1) = L_pieces.*cast( tf, 'like', L.pieces );
%                         figure(1); subplot( 1, 2, 1 ); cla; imagesc( L_conn ); title( num2str( p' ) )
                    
                    bw_conn = comp_reg( L_conn( [1 1:end-1], : ), L_conn( [2:end end], : ) ) | ...
                        comp_reg( L_conn( :, [1 1:end-1] ), L_conn( :, [2:end end] ) ) | ...
                        comp_reg( L_conn( [1 1:end-1], [1 1:end-1] ), L_conn( [2:end end], [2:end end] ) ) | ...
                        comp_reg( L_conn( [2:end end], [1 1:end-1] ), L_conn( [1 1:end-1], [2:end end] ) );
                    bw_conn = bw_conn(2:end-1, 2:end-1)|tf;
                    
                    %%
                    %%
                    if sum( bw_conn(:) ) >= cutoff.Area
                        CC.pieces.NumObjects = CC.pieces.NumObjects + 1;
                        [r, c] = find( bw_conn );
% %                         figure(1); subplot( 1, 2, 1 ); cla; imagesc( bw_conn );
                        if numel( CC.pieces.PixelIdxList ) < CC.pieces.NumObjects
                            CC.pieces.PixelIdxList{10*end} = [];
                        end
                        CC.pieces.PixelIdxList{CC.pieces.NumObjects} = sub2ind( CC.input.ImageSize, r+r_lim(1)-1, c+c_lim(1) - 1 );
% %                        figure(1); subplot( 1, 2, 2 ); cla;  tmp = false( CC.pieces.ImageSize ); tmp( CC.pieces.PixelIdxList{CC.pieces.NumObjects} ) = true; imagesc( tmp )
                        if numel( CC.pieces.MaxCutLength ) < CC.pieces.NumObjects
                            CC.pieces.MaxCutLength( 10*end) = 0;
                            CC.pieces.CutLength( 10*end) = 0;
                            CC.pieces.level( 10*end) = 0;
                        end
                        a = unique( vertcat( Neighbors{p} ) );
                        b = setdiff( a, p );
                        if ~isempty( b )
                            CC.pieces.MaxCutLength(CC.pieces.NumObjects) = max( reshape( cut_length_pieces12( a, b ), [], 1 ) );
                            CC.pieces.CutLength(CC.pieces.NumObjects) = sum( reshape( cut_length_pieces12( a, b ), [], 1 ) );
                        end
                        [CC.pieces.level(CC.pieces.NumObjects), level(CC.pieces.PixelIdxList{CC.pieces.NumObjects})] = ...
                            deal( max( level(CC.pieces.PixelIdxList{CC.pieces.NumObjects}) ) + 1 );
                        %                         imagesc( level ); hold on
                    end
                    %                             cla; imagesc( ( L.input > 0 ) + (L.pieces > 0 ) ); hold on; [y, x] = ind2sub( CC.input.ImageSize, CC.pieces.PixelIdxList{CC.pieces.NumObjects} ); plot( x, y, '.r' ); title( num2str( p ) )
                elseif false
                    L_conn = L.pieces.*cast( ismember( L.pieces, p ), 'like', L.pieces );
                    bw_conn = comp_reg( ...
                        L_conn( [1 1:end-1], : ), L_conn( [2:end end], : ) ) | comp_reg( L_conn( :, [1 1:end-1] ), L_conn( :, [2:end end] ) ) | ...
                        comp_reg( L_conn( [1 1:end-1], [1 1:end-1] ), L_conn( [2:end end], [2:end end] ) ) | ...
                        comp_reg( L_conn( [2:end end], [1 1:end-1] ), L_conn( [1 1:end-1], [2:end end] ) ) | ...
                        ismember( L.pieces, p );
                    cla; imagesc( ( L.input > 0 ) + (L.pieces > 0 ) ); hold on; [y, x] = find( bw_conn ); plot( x, y, '.r' ); title( num2str( p, '*%g ' ) )
                end
            end
        end
    end
end
%%
[CC.pieces.level, CC.pieces.MaxCutLength, CC.pieces.CutLength, CC.pieces.PixelIdxList] = ...
    deal( CC.pieces.level(1:CC.pieces.NumObjects ), CC.pieces.MaxCutLength(1:CC.pieces.NumObjects ), CC.pieces.CutLength(1:CC.pieces.NumObjects ), CC.pieces.PixelIdxList(1:CC.pieces.NumObjects ) );
% [CC.pieces.Area, ind] = sort( cellfun( @numel, CC.pieces.PixelIdxList ), 'descend' );
% CC.pieces.PixelIdxList = CC.pieces.PixelIdxList(ind);
varargout = {CC.pieces};
% [Area, ind] = sort( Area, 'descend' );
% [IdxFill, IdxEdge, Mother, Daughter] = deal( IdxFill(ind), IdxEdge(ind), Mother(ind), Daughter(ind) );
%%
if nargout == 0
    fig = figure( sum( mfilename ) );
    set( fig, 'Name', mfilename, 'UserData', struct( 'CC', CC.pieces, 'level', 1, 'max_level', max( CC.pieces.level ) ), 'KeyPressFcn', @KeyPressFcn );
    sp = [subplot( 1, 2, 1, 'NextPlot', 'add', 'DataAspectRatio', [1 1 1], 'Parent', fig, 'Tag', 'input' ), subplot( 1, 2, 2, 'NextPlot', 'add', 'DataAspectRatio', [1 1 1], 'Parent', fig, 'Tag', 'pieces' )];
    cla( sp(1) );
    imagesc( sp(1), labelmatrix( CC.input ) );
    linkaxes( sp );
    axis tight;
    update_plot( fig )
end

function KeyPressFcn( fig, event )

switch event.Key
    case 'uparrow'
        fig.UserData.level = fig.UserData.max_level;
    case 'downarrow'
        fig.UserData.level = 1;
    case 'leftarrow'
        fig.UserData.level = mod( fig.UserData.level-2, fig.UserData.max_level ) + 1;
    case 'rightarrow'
        fig.UserData.level = mod( fig.UserData.level, fig.UserData.max_level ) + 1;
end
update_plot( fig );

function update_plot( fig )
sub = find( fig.UserData.CC.level == fig.UserData.level );
ax = findobj( fig, 'Tag', 'pieces' );
img = findobj( ax, 'Type', 'image' );
if isempty( img )
    imagesc( ax, labelmatrix_sak( fig.UserData.CC, sub ), 'ButtonDownFcn', @ButtonDownFcn );
else
    set( img, 'CData', labelmatrix_sak( fig.UserData.CC, sub ) );
end
title( ax, sprintf( 'level %i', fig.UserData.level ) );
%%

function bnry = get_binary( n )
persistent cache
if isempty( cache ) ||  n > size( cache, 2 )
    cache = fliplr( dec2bin( 1:2^n-1 ) == '1' );
end
bnry = cache( 1:2^n-2, 1:n );

function ButtonDownFcn( img, event )
%%
ind = img.CData( round( event.IntersectionPoint(2) ), round( event.IntersectionPoint(1) ) );
if ind > 0
    fprintf( '# %i, CutLength: %i, MaxCutLength: %i\n', ind, img.Parent.Parent.UserData.CC.CutLength( ind ),img.Parent.Parent.UserData.CC.MaxCutLength( ind ) );
end
