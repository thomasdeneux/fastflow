
% what to do
createraw = 1;
dopreproc = 1;
doanalysis = 0;

% main directory where all the data and results will be saved
root_dir = '/sulci3/deneux/Projects/0503_Marseille/BiasProblem/fakedata64normalizeagain';

% list of speeds (of the RBC) to be tested
step=1.5; %[1 .5 2];

% list of thresholds to be tested, which controls the density
% of the RBCs: the higher the threshold, the less RBCs
thresh=[0.8];

% list of values for the jitter
jitter=[0.6 0];

% list of values for the binning (1 = no binning)
binning=[1];

% number of simulated cases for one set of parameters defined above
simul_nbr_per_case = 64;

% parameters of the simulated datatsets
nj=100; % y-dimension: the RBCs shall move in a vessel oriented vertically
ni=100; % x-dimension
nt=300; % number of frames

% parameters for the log-gabor filter set
%nbr_scales = 4;
nbr_orient = 89;
min_angle = 1; % in degrees
max_angle = 89; % in degrees
angles = min_angle + (0:nbr_orient-1)/(nbr_orient-1)*(max_angle-min_angle);
factper = power(1.3,[0 1 2]);
factext = power(2,[0 1 2]);
periods = 4 * factper([1 1 1 2 2 2 3 3 3]);
extents = periods .* (1/2 * factext([1 2 3 1 2 3 1 2 3]));
ellips = 1/8 * extents;
filter_size = pow2(floor(log2(max([ni nj nt]))) + 1);
% clear gaborCoef

% compute gabor filter that will be used to determine in which quarter we're working!
% use only 2 orientations (-45 degrees and +45), and 3 scales to make it robust
quarter_choice_filter_struct = log_gabor_pyramid_wedge(filter_size,3,2,-45*pi/180,45*pi/180);

% create raw datasets
load diversTOTALfake
load paramsfake
edges=edges(1);


for s=1:length(step)  % loop over the speed values
    for t=1:length(thresh)  % loop over the threshol values
        % go to the good directory and create subdir
        cd(root_dir);
        sub_dir1 = sprintf('step_%03d_thresh_%03d',100*step(s),100*thresh(t));
        if ( exist(sub_dir1)==0 )
            mkdir(sub_dir1);
        end;
        full_sub_dir1 = fullfile(root_dir,sub_dir1,'');
        
        for k = 1:simul_nbr_per_case  % loop of simuls
            
            fprintf('\nSPEED %g\nTHRESHOLD %g\nSIMUL %i\n',step(s),thresh(t),k)
            % CREATE RAW DATA
            if createraw
                disp('Create raw data'), tic
                % create raw data (this function now acts once for all the
                % values of jitters specified
                fake_all = NEWFAKE_syl(thresh(t),step(s),jitter,ni,nj,nt);
                toc
            end
            
            for j=1:length(jitter)  % loop over jitter values
                for b=1:length(binning)  % loop over binning values
                    
                    fprintf('JITTER %g\nBINNING %g\n',jitter(j),binning(b))
            
                    % create new subdir and get in there
                    cd(full_sub_dir1)
                    sub_dir2 = sprintf('jitter_%03d_binning_%02d',100*jitter(j),binning(b));
                    if ( exist(sub_dir2,'dir')==0 )
                        mkdir(sub_dir2);
                    end;
                    full_sub_dir2 = fullfile(full_sub_dir1,sub_dir2,'');
                    cd(full_sub_dir2);
                    
                    % SAVE RAW DATA
                    if createraw
                        fake=fake_all{j};
                        simul_name = sprintf('simul_%03d.mat',k);
                        save(fullfile(full_sub_dir2,simul_name),'fake');
                    end
                    
                    % PREPROC
                    if dopreproc
                        if ~createraw
                            simul_name = sprintf('simul_%03d.mat',k);
                            load(fullfile(full_sub_dir2,simul_name),'fake')
                        end
                        
                        disp('Preprocessing'), tic
                        CS=fake.CS;
                        CSU=fake.CSU;
                        CSV=fake.CSV;
                        CS0 = CS - mean(CS(:));
                        Y_mean = fake.Y;
                        %Y_mean = fast_normalizeseq1test_mean('Y_mean',binning(b)); %original by SylvainT
                        
                        %Y_mean = fast_normalizeseq('Y_mean',binning(b)); 
                        Y_mean = fast_normalize('Y_mean');
                        
                        %Y = fast_highpass('Y',10); disp(datestr(now,'HH:MM:SS'))
                        
                        %Y = fast_anisosmooth('Y',CSU(7:344,7:220),CSV(7:344,7:220)); 
                        %disp(datestr(now,'HH:MM:SS'))
                        
                        
                        % another binning method!
                        % I havn't tested what the influence is on the angle
                        % estimation yet... you should do it, ivo!
                        %Y_min  = fake.Y;
                        %Y_min  = fast_normalizeseq1test_min('Y_min',binning(b)); 
                        
                        
                        preproc_name = sprintf('preproc_%03d.ff',k);
                        fast_save2bytes(Y_mean,preproc_name)
                        toc
                    end
                    
                    % ANALYSIS
                    if doanalysis
                        if ~dopreproc
                            preproc_name = sprintf('preproc_%03d.ff',k);
                            Y_mean = fast_load2bytes(preproc_name);
                        end
                        
                        disp('Analysis'), tic
                        e = fast_interpalongedges(edges,Y_mean);

                        if ~exist('gaborCoef','var')
                            % estimate in which quarter we're working; returns +1
                            % or -1
                            quarter = select_quarter(e,quarter_choice_filter_struct);
                            if quarter == 0
                                error(['Problem: I cannot determine in which ' ...
                                        'quarter we are working']);
                            end;
                            % adds 90 degres to the limit angles for the
                            % given quarter
                            if quarter == -1
                                min_angle = min_angle + 90;
                                max_angle = max_angle + 90;
                            end;
                            % now, define the full gabor analysis filter
                            % this assumes all the simulated data "live" in
                            % the same quarter
                            gaborCoef = TOM_gaborbase(filter_size,periods,extents,ellips,angles*pi/180);
                        end;
                        disp('apply Gabor')
                        res = TOM_applygabor(e,1,0,gaborCoef);
                        res_name =  sprintf('res_%03d.mat',k);
                        save(res_name,'res');
                        toc
                    end
                    
                end;
            end;
        end;
    end;
end;

