function f = get_hem_folder()
%HEM.UTIL.GET_HEM_FOLDER Get +hem folder
%   F = HEM.UTIL.GET_HEM_FOLDER returns the full path to the +hem folder.
%   This is useful because HEM.MODEL.RUN (and others) need to know the
%   location of the +hem/params folder to load parameter values.
%
%   If there are multiple +hem folders MATLAB's current search path, the
%   first one is selected. So be careful!

w = what('hem');
f = w(1).path;

if length(w) > 1
    wstr = 'Multiple +hem folders found in MATLAB''s current search path:\n';
    for i=1:length(w)
        wstr = [wstr, sprintf('%s\n', w(i).path)];
    end
    wstr = [wstr, sprintf('Currently using the first path: %s', f)];
    warning('hem:multiple_folders',wstr)
end
