% Raw2Stl.m
%
% Babak Taati 
% http://qlink.queensu.ca/~3bt1/
% Queen's University
% Feb 16, 2005
%
% a function for converting 3d objects of RAW format to STL format 
%
% syntax:
%   boolean b = Raw2Stl(strRawFileName, strStlFileName);
%   
% returns 'false' when fails 
%
% possible failture reasons:
%   - input file is not of raw format
%   - output file cannot be created
%
% Note: normal sign ambiguity ... (information in RAW file not enough to resolve the sign ambiguity)
%

function b = Raw2Stl(RawData, strStlFileName)

if ischar(RawData)
    RawData = load(RawData); % should be an  n*9  matrix
end
    
% note: the 'load' function fails if the file cannot be read (e.g. if the number of columns is not the
%       same in all rows)


[NoOfTriangles, ShouldBeNine] = size(RawData);

if(  (ShouldBeNine ~= 9)  ||  (NoOfTriangles == 0)   )
    fprintf('Error: input file not in RAW format!');
    b = false;
    return;
end

P1 = RawData(:,1:3); % P1(ii,:) will be the 1st vertex of the ii'th triangle
P2 = RawData(:,4:6); % P2(ii,:) will be the 2nd vertex of the ii'th triangle
P3 = RawData(:,7:9); % P3(ii,:) will be the 3rd vertex of the ii'th triangle

U = P2 - P1;    % a side of each triangle
V = P3 - P1;    % the other side of each triangle

TriangleNormals = cross(U,V);   %   1- sign ambiguity ... (can't be resolved)
                                %   2- yet to be normalized
                                

% --start: write to file

fid = fopen(strStlFileName, 'w');
                                
if(fid == -1)
    fprintf('Error: could not write to file');
    b = false;
    return;    
end

fprintf(fid, 'solid Object\n');

for ii = 1:NoOfTriangles
    
    NormalizedTriangleNormal = TriangleNormals(ii,:) / norm( TriangleNormals(ii,:) );

    fprintf(fid, ' facet normal %f %f %f\n', NormalizedTriangleNormal);

    fprintf(fid, '  outer loop\n');

    fprintf(fid, '   vertex %f %f %f\n', P1(ii,:) );
    fprintf(fid, '   vertex %f %f %f\n', P2(ii,:) );
    fprintf(fid, '   vertex %f %f %f\n', P3(ii,:) );    

    fprintf(fid, '  endloop\n');
    
    fprintf(fid, ' endfacet\n');
   
end

fprintf(fid, 'endsolid\n');

fclose(fid);

% --end: write to file

b = true;
return