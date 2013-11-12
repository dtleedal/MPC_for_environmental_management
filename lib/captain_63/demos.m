function tbxStruct=demos
% DEMOS  Returns demo information

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

if nargout==0, demo toolbox; return; end

tbxStruct.Name='Captain';
tbxStruct.Type='toolbox';

tbxStruct.Help={
         ' The Captain Toolbox is a collection of         '  
         ' functions for non-stationary time series       '
         ' analysis and forecasting.                      '
         '                                                '
         ' It allows for the identification of unobserved '
         ' components models, time variable parameter     '
         ' models, state dependent parameter models and   '
         ' multiple-input discrete and continuous time    '
         ' transfer function models.                      '
         '                                                '
         ' The toolbox is useful for system identification, '
         ' signal extraction, interpolation, backcasting, '
         ' forecasting and data-based mechanistic analysis'
         ' of a wide range of linear and non-linear       '
         ' stochastic systems.'};

tbxStruct.DemoList={
         'Background', 'captinfo',
         'Air Passengers', 'playshow dhrshow',
         'Chirp Signal', 'playshow darshow',
         'Transfer Function', 'playshow rivshow',
         'Command Line Demos 1', 'captdems1',
         'Command Line Demos 2', 'captdems2',
         'Command Line Demos 3', 'captdems3'};
