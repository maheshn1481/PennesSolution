% MPH2STL Convert COMSOL geometry object to stl-file
%   requires Raw2stl.m by Babak Taati (slightly modified to accept matrix
%   instead of file)
%     MPH2STL(infile,outfile) loads the mphbin/mphtxt-file, plots it and
%     calls Raw2stl to save the stl file. Plotting is essential because the
%     triangulated geometry is taken from the plot
%     MPH2STL(infile) does the same as above, but saves the data in
%     (infile without extension)+'stl'
%     MPH2STL or MPH2STL([],outfile) searches the workspace for a geometry
%     object (exported from COMSOL in linked mode) and converts it
% 
%   Henning Francke 2007-04-05.
%   Copyright (c) 2007 by Weierstraﬂ Institute Berlin (www.wias-berlin.de)
function mph2stl(infile,outfile)

if ~exist('infile','var')||isempty(infile)
    if ~exist('outfile','var')
        outfile = 'mph2stl.stl';
    end
    disp('No filename passed, searching for COMSOL geometry object in workspace...')
    ws = evalin('base','whos');
    for i=1:length(ws)
        if strcmp(ws(i).class,'solid3')
            disp('...found!')
            gemalt = geomplot(evalin('base',ws(i).name));
    %             o = eval(ws(i).name);
            break
        end
    end
    if ~exist('gemalt','var')
        disp('No geometry object found in workspace!')
        return
    end
else
    objekte = geomimport(infile);
    gemalt = geomplot(objekte{1});
    if ~exist('outfile','var')
        outfile = [regexprep(infile,'.mph(bin|txt)','') '.stl'];
    end
end

% kinder = get(achse,'children');
% kinder = get(gca,'children');
% v = flgeomvtx(o);    %Get vertices
% nf = flgeomnbs(o);
% [xyz,dxyz]=flgeomfd(EXT1,1:nf,repmat([0;1],1,nf));
% n = flgeomfn(dxyz);  %Normale

% m = com.femlab.mesh.GeomMeshUtil.geomMeshML(getjptr(o),'normal',0,1);
% mesh = flgeommesh(o,'normal',0,1);
% p=mesh.p;
% eg=mesh.eg;
% e=mesh.e;
% t=mesh.t;

% fid = fopen('test.stl', 'w');
% fprintf(fid,'solid Object\n');
rawdata=[];
for i=1:length(gemalt)
    if strcmp(get(gemalt(i),'type'),'patch')        
        v = get(gemalt(i),'vertices');
        f = get(gemalt(i),'faces');
        for j=1:size(f,1)
%             fprintf(fid,' facet normal 0.000000 0.000000 1.000000\n');
%             fprintf(fid,'  outer loop\n');
%                 for k=1:size(f,2)
%                      fprintf(fid,'   vertex %.5f %.5f %.5f\n',v(f(j,k),:));
%                 end
                    rawdata(end+1,:)= reshape(v(f(j,:),:)',1,[]);
%            fprintf(fid,'  endloop\n endfacet\n');
        end
    end
end
if Raw2stl(rawdata,outfile)
    disp(['Geometry apparently successfully converted and saved to "' outfile '".'])
end
% fprintf(fid,'endsolid');
% fclose(fid);