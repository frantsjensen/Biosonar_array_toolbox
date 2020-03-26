function [group] = separateclickseries(C)
% function [group] = separateclickseries(C)
% 
% Separate overlapping odontocete click trains based on ICI
% C is a vector of time stamps for each click (in samples or seconds)
% Group is a vector of group identity for each click
%
% Separateclicks is not very reliable and should be used in conjunction
% with visualization. Currently supports two possible groups.
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

% Default group for each click
group=ones(1,length(C));

% Find clicks that are in doubt
kdoubt = find(diff(C)<0.5*median(diff(C)));

% If no clicks below ICI threshold, break script
if isempty(kdoubt)
    return
end

% Otherwise find potential trouble clicks
kdoubt = [ kdoubt(:) kdoubt(:)+1] ;
kgood  = setxor([1:length(C)],kdoubt(:)) ;

% Construct matrix with new index, old index, and observed time for each good click
M = sortrows([kdoubt(:,1) zeros(size(kdoubt,1),1) ; kgood(:) C(kgood)']);
M = [ [1:length(M)]' M];

% Construct matrix with proximity to most observations, and old index, and
% sort after proximity
prox = sortrows([abs(C(kdoubt(:,1))'-mean(C(kgood))) , [1:size(kdoubt,1)]' ]);

% Now work through clicks that are in doubt, starting from click closest to
% most of our known points
for j=1:size(prox,1)
    oldindex = kdoubt(prox(j,2)) ;
    krow = find(M(:,2)==oldindex) ;
    kgood = find(M(:,3)>0) ;
    Cobs = C(kdoubt(prox(j,2),:));
    Cest = interp1(M(kgood,1),M(kgood,3),krow,'spline','extrap');
    [test,correct] = min(abs(Cobs-Cest)) ;
    [test,wrong] = max(abs(Cobs-Cest)) ;
    k2(j) = kdoubt(prox(j,2),wrong) ;
    M(krow,3) = C(kdoubt(prox(j,2),correct)) ; % Update matrix with known click sample
end

group(sort(k2)) = 2;


% ALTERNATIVES (NOT VERY GOOD)

% % Always start with group 1
% currentgroup = 1;
% group=zeros(1,length(C));
% 
% % Iteratively go through click series to search for possible mixed scans
% for j=1:length(C)
%     group(j) = currentgroup ; % Set group for this click
%     switch j
%         case 1
%             ind  = [1 2 3] ; 
%             ind2 = [1 3 5] ;
%         case 2
%             ind  = [1 2 3] ; 
%             ind2 = [2 4 6] ;
%         case length(C)-1
%             ind  = [-2 -1 0]+length(C)-1 ;
%             ind2 = [-4 -2 0]+length(C)-1 ;
%         case length(C)
%             ind  = [-2 -1 0]+length(C) ;
%             ind2 = [-4 -2 0]+length(C) ;
% 
%         otherwise
%             ind  = [-1 0 1]+max([3 min([j length(C)-2]) ]) ; 
%             ind2 = [-2 0 2]+max([3 min([j length(C)-2]) ]) ;
%     end
%     
%     % Then see if next click should be in different group
%     if std(diff(C(ind2))) < 0.5*std(diff(C(ind))),
%         currentgroup = 3-currentgroup ; % Alternate between two groups
%     end
% end



% % Starkhammar suggested correlating power spectra
% for j=1:length(Ct),
%     thissig = sig(j,:) .* tukeywin(length(sig(j,:)),0.15)';
%     sigfft  = abs(fft(thissig,256)) ;
%     sigfft  = sigfft(1:length(sigfft)/2) ;
%     sigfft  = sigfft/max(sigfft);
%     specpow(j,:) = sigfft ;
% end
%
% % But this does not work
% thesefreqs = 29:length(specpow);
% for j=1:length(Ct)-2,
%     CC(j) = sum( specpow(j,thesefreqs).*specpow(j+1,thesefreqs) )/ max([sum(specpow(j,thesefreqs).^2) sum(specpow(j+1,thesefreqs).^2)]) ;
% end
