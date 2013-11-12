function forcing = IAGPgetForcing(value,species,units)%CO2_in,CH4_in,N2O_in,CFC_11_in,CFC_12_in)
% forcing = getForcing(value,species,units)
% forcing:  the GHG forcing for the given input (W m^-2)
% value:    the concentration of the input gas
% species:  a string describing the input gas, must be one of
%           ['CO2'], 'CH4', 'N2O', 'CFC11', 'CFC12', 'CFC13', 'HFC134A','O3'
% units:    the units of the input gas concentration, must be one of
%           ['ppm'], 'ppb'

if nargin == 1
    species = 'CO2';
    disp('Assuming input is CO2')
    units = 'ppm';
    disp('Assuming units is ppm')
end
if nargin == 2
    units = 'ppm';
    disp('Assuming units is ppm')
end

if isempty(find(strcmp(species,{'CO2', 'CH4', 'N2O', 'CFC11', 'CFC12', 'CFC13', 'HFC134A','O3'}))==1)
error('string describing the input type must be one of: ''CO2'', ''CH4'', ''N2O'', ''CFC11'', ''CFC12'', ''CFC13'', ''HFC134A'',''O3''')
end

if isempty(find(strcmp(units,{'ppm', 'ppb'}))==1)
error('string describing the input units must be one of: ''ppm'', ''ppb''')
end

% pre industrial value
CO2_0 = 278.6;% in ppm
CH4_0 = 722;% in ppb
N2O_0 = 273;% in ppb
CFC_11_0 = 0;% in ppb
CFC_12_0 = 0;% in ppb
CFC_13_0 = 0;% in ppb

% perform the various conversions
% ---------------------CO2-----------------------
if strcmp(species,'CO2')
    % do CO2 conversion
    if strcmp(units,'ppb')
        % convert to ppm
        CO2 = value * 1e-3;
    else
        CO2 = value;
    end
    % pre industrial values (_0 appended)
    
    %CO2 forcing
    forcing = 5.35.*log(CO2./CO2_0);
end
% -------------------Methane---------------------
if strcmp(species,'CH4')
    % do methane conversion
    if strcmp(units,'ppm')
        % convert to ppb
        CH4 = value * 1e3;
    else
        CH4 = value;
    end      
        %Methane forcing
% forcing = 0.036*(sqrt(CH4)-sqrt(CH4_0)) ...
%     - 0.47*log(1 + 2.01e-5*((CH4*N2O_0).^0.75) +... 
% 5.31e-15*CH4.*((CH4*N2O_0).^1.52)) ...
%     + 0.47*log(1 + 2.01e-5*((CH4_0*N2O_0).^0.75) + 5.31e-15*CH4_0.*((CH4_0*N2O_0).^1.52));
forcing = 0.036.*(sqrt(CH4) - sqrt(CH4_0)) - ( f(CH4,N2O_0) - f(CH4_0,N2O_0) );

end
% --------------------NOX------------------------
if strcmp(species,'N2O')
    if strcmp(units,'ppm')
        % convert to ppb
        N2O = value .* 1e3;
    else
        N2O = value;
    end
    % pre industrial value
    
%     forcing = 0.12*(sqrt(N2O)-sqrt(N2O_0)) ...
%     - 0.47*log(1 + 2.01e-5*((CH4_0*N2O).^0.75) +... 
% 5.31e-15*CH4_0.*((CH4_0*N2O).^1.52)) ...
%     + 0.47*log(1 + 2.01e-5*((CH4_0*N2O_0).^0.75) + 5.31e-15*CH4_0.*((CH4_0*N2O_0).^1.52));
forcing = 0.12.*(sqrt(N2O) - sqrt(N2O_0)) - ( f(CH4_0,N2O) - f(CH4_0,N2O_0) );
end
% --------------------CFC11----------------------
if strcmp(species,'CFC11')
    if strcmp(units,'ppm')
        % convert to ppb
        CFC_11 = value .* 1e3;
    else
        CFC_11 = value;
    end
    % pre industrial
    
    forcing = 0.25.*(CFC_11 - CFC_11_0);
end
% --------------------CFC12----------------------
if strcmp(species,'CFC12')
    if strcmp(units,'ppm')
        % convert to ppb
        CFC_12 = value .* 1e3;
    else
        CFC_12 = value;
    end
    % pre industrial
    
    forcing = 0.32.*(CFC_12 - CFC_12_0);
end
% --------------------CFC13----------------------
if strcmp(species,'CFC12')
    if strcmp(units,'ppm')
        % convert to ppb
        CFC_13 = value * 1e3;
    else
        CFC_13 = value;
    end
    % pre industrial
    
    forcing = 0.25.*(CFC_13 - CFC_13_0);
end
% ----------------------O3-----------------------

if strcmp(species,'O3')
    error('Not able to process O3 values at present WTS')
end
end

function val = f(M,N)
val = 0.47 .* log( 1 + 2.01e-5 .* (M.*N).^0.75 + 5.31e-15 .* M.*(M.*N).^1.52 );
end




