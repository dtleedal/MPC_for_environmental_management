% script to get input and run IAGPperformMPC


disp('                         * * * * * * * ')
disp('                         *  WELCOME  *')
disp('                         * * * * * * *')
disp('  ')
disp('Step 1, confirm you are in the MPC_for_Lawrence directory')


this_directory = pwd;
disp('The current directory is:')
disp('  ')
disp(pwd)
disp('  ')
answer = input('Make sure this is the right directory. If it is then enter [y]. Otherwise enter \nany other letter to quit, then change to the right directory and run this file again\n>','s');
if answer~='y'
    return
end
confirm_response = 'n';
while confirm_response ~='y'
    disp('Step 2, Collecting the input data')
    
    disp('  ')
    year = input('Enter the year that has just finished: ');
    NH = input('Enter the observed Northern hemisphere temperature for that year: ');
    SH = input('Enter the observed Southern hemisphere temperature for that year: ');
    global_temp = input('Enter the observed global temperature for that year: ');
    sea_ice_min = input('Enter the observed minimum sea ice extent for that year: ');
    emis = input('Enter the SO2 emissions (Tg) that were applied in that year: ');
    ghg_forcing = input('Enter the GHG forcing that was applied in that year (not including volcano effects): ');
    volcano_size = input('Enter the approximate SO2 injection (Tg) if there was a volcanic event last year: ');
    disp('  ')
    disp('Check your input:')
    disp(['Year:               ' num2str(year)])
    disp(['NH temperature:     ' num2str(NH)])
    disp(['SH temperature:     ' num2str(SH)])
    disp(['global temperature: ' num2str(global_temp)])
    disp(['sea ice minimum:    ' num2str(sea_ice_min)])
    disp(['SO2 emissions:      ' num2str(emis)])
    disp(['the GHG forcing:    ' num2str(ghg_forcing)])
    disp(['volcanic event:     ' num2str(volcano_size)])
    disp('   ')
    if volcano_size > 0
        disp({('There has been a volcano! continue on with this script');...
            ('but then request a re-spin before running HadGEM')});
    end
    confirm_response = input('Is this correct? Enter [y] if it is: ','s');
end

disp('  ')
confirm_response = 'blank';
while confirm_response~='r' | confirm_response~='a'
    confirm_response = input('Step 3, hit [r] to run [b] to bail out: ','s');
    
    if confirm_response == 'r'
        disp('Running...')
        emissions = IAGPperformMPC(year,[NH; SH; global_temp; sea_ice_min],...
            emis,ghg_forcing,volcano_size,this_directory);
        return
    end
    if confirm_response == 'b'
        disp('...bailing out!')
        return
    end
end


