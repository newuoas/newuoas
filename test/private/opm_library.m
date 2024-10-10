function opm_library()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Describes and verifies all problems in the OPM library.
%
%   Programming : Ph. Toint, July 2018
%   This version: 28 VII 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oneline = 1;

disp( ' ' )
disp(   '*************************************************************************************************************************'  )
disp(   '*************************************************************************************************************************'  )
disp(   '********************************                                           **********************************************'  )
disp( [ '********************************     OPM library content (', date, ')     **********************************************' ] )
disp(   '********************************                                           **********************************************'  )
disp(   '*************************************************************************************************************************'  )
disp(   '*************************************************************************************************************************'  )
disp( ' ' )
fprintf( 'Problem name        fstar         n    nc   ni nfree  nlow  nupp nboth  nfix  m   mi   me objtype constype  classif\n' )
disp( ' ' )

opm_describe( 'argauss'      , oneline );
opm_describe( 'arglina'      , oneline );
opm_describe( 'arglinb'      , oneline );
opm_describe( 'arglinc'      , oneline );
opm_describe( 'argtrig'      , oneline );
opm_describe( 'arwhead'      , oneline );
opm_describe( 'bard'         , oneline );
opm_describe( 'bdarwhd'      , oneline );
opm_describe( 'beale'        , oneline );
opm_describe( 'biggs5'       , oneline );
opm_describe( 'biggs6'       , oneline );
opm_describe( 'booth'        , oneline );
opm_describe( 'box2'         , oneline );
opm_describe( 'box3'         , oneline );
opm_describe( 'bratu1d'      , oneline );
opm_describe( 'brkmcc'       , oneline );
opm_describe( 'brownbs'      , oneline );
opm_describe( 'brownden'     , oneline );
opm_describe( 'brownal'      , oneline );
opm_describe( 'broyden3d'    , oneline );
opm_describe( 'broydenbd'    , oneline );
opm_describe( 'chebyqad'     , oneline );
opm_describe( 'cliff'        , oneline );
opm_describe( 'clustr'       , oneline );
opm_describe( 'cosine'       , oneline ); %
opm_describe( 'crglvy'       , oneline );
opm_describe( 'cube'         , oneline ); %
opm_describe( 'dixmaana'     , oneline );
opm_describe( 'dixmaanb'     , oneline );
opm_describe( 'dixmaanc'     , oneline );
opm_describe( 'dixmaand'     , oneline );
opm_describe( 'dixmaane'     , oneline );
opm_describe( 'dixmaanf'     , oneline );
opm_describe( 'dixmaang'     , oneline );
opm_describe( 'dixmaanh'     , oneline );
opm_describe( 'dixmaani'     , oneline );
opm_describe( 'dixmaanj'     , oneline );
opm_describe( 'dixmaank'     , oneline );
opm_describe( 'dixmaanl'     , oneline );
opm_describe( 'dixon'        , oneline );
opm_describe( 'dqrtic'       , oneline );
opm_describe( 'edensch'      , oneline );
opm_describe( 'eg2'          , oneline);
opm_describe( 'eg2s'         , oneline );
opm_describe( 'engval1'      , oneline );
opm_describe( 'engval2'      , oneline );
opm_describe( 'freuroth'     , oneline );
opm_describe( 'genhumps'     , oneline );
opm_describe( 'gottfr'       , oneline );
opm_describe( 'gulf'         , oneline );
opm_describe( 'hairy'        , oneline ); %
opm_describe( 'heart6ls'     , oneline ); %
opm_describe( 'helix'        , oneline );
opm_describe( 'hilbert'      , oneline );
opm_describe( 'himln3'       , oneline );
opm_describe( 'himm25'       , oneline );
opm_describe( 'himm27'       , oneline );
opm_describe( 'himm28'       , oneline );
opm_describe( 'himm29'       , oneline );
opm_describe( 'himm30'       , oneline );
opm_describe( 'himm32'       , oneline );
opm_describe( 'himm33'       , oneline );
opm_describe( 'hypcir'       , oneline );
opm_describe( 'indef'        , oneline );
opm_describe( 'integreq'     , oneline ); %
opm_describe( 'jensmp'       , oneline );
opm_describe( 'kowosb'       , oneline );
opm_describe( 'lminsurf'     , oneline );
opm_describe( 'mancino'      , oneline ); %
opm_describe( 'mexhat'       , oneline ); %
opm_describe( 'meyer3'       , oneline ); %
opm_describe( 'morebv'       , oneline );
opm_describe( 'msqrtals'     , oneline );
opm_describe( 'msqrtbls'     , oneline );
opm_describe( 'nlminsurf'    , oneline );
opm_describe( 'osbornea'     , oneline ); %
opm_describe( 'osborneb'     , oneline ); %
opm_describe( 'penalty1'     , oneline );
opm_describe( 'penalty2'     , oneline );
opm_describe( 'penalty3'     , oneline ); %
opm_describe( 'powellbs'     , oneline ); %
opm_describe( 'powellsg'     , oneline ); %
opm_describe( 'powellsq'     , oneline ); %
opm_describe( 'powr'         , oneline );
opm_describe( 'rosenbr'      , oneline ); %
opm_describe( 'schmvett'     , oneline );
opm_describe( 'scosine'      , oneline );
opm_describe( 'sisser'       , oneline );
opm_describe( 'spmsqrt'      , oneline ); %
opm_describe( 'recipe'       , oneline );
opm_describe( 'tridia'       , oneline );
opm_describe( 'trigger'      , oneline );
opm_describe( 'vardim'       , oneline );
opm_describe( 'watson'       , oneline );
opm_describe( 'wmsqrtals'    , oneline );
opm_describe( 'wmsqrtbls'    , oneline );
opm_describe( 'woods'        , oneline );
opm_describe( 'zangwil2'     , oneline );
opm_describe( 'zangwil3'     , oneline );

disp( ' ' )
disp(   '*************************************************************************************************************************'  )
disp( ' ' )
