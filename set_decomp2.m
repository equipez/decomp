function [P2, Pc2, mask2, ratio_ff] = set_decomp2 (n, m, n_o, as_type, coarse, stride)
% Gratton Zhang 2017
% 1D input, but 2D components output
% n number of variables
% m number of subspaces
% n_o size of the overlap
% as_type 'as' or 'ras' or ?? for weighting
% each is 1D. Call set_decomp to build a 2D decomp.
 
%[P0, Pc0, mask0] = set_decomp (n, m, 0, as_type);
%[P0, Pc0, mask0] = set_decomp (n, m, n_o, as_type);
%[P, Pc, mask] = set_decomp (n, m, n_o, as_type);
P0_1 = set_decomp (n, m(1), n_o, as_type);
P0_2 = set_decomp (n, m(2), n_o, as_type);
[P_1, ~, mask_1] = set_decomp (n, m(1), n_o, as_type);
[P_2, ~, mask_2] = set_decomp (n, m(2), n_o, as_type);

ratio_ff=1;

index2d =reshape(1:n^2,n,n);
II = speye(n^2,n^2);
%
ind2 = 1;
last_dom={};
Pc2=[];
for i=1:m(1)
  for j=1:m(2)
    coarse_handeled = sparse(zeros(n,n));% the coarse point has been analyzed
    indP2{ind2} = nonzeros(index2d.*(sum(P_1{i},2)*sum(P_2{j},2)')) ; % list of all nonzeros entries
    P2{ind2} = II(:,indP2{ind2});
    mask2{ind2}=reshape(mask_1{i}*mask_2{j}',1,length(mask_1{i})*length(mask_2{j}'))';
if (coarse == 6 | coarse == 8)
% Coarse grid fonstruction
xmin=min(find(sum(P0_1{i},2)));xmax=max(find(sum(P0_1{i},2)));
ymin=min(find(sum(P0_2{j},2)));ymax=max(find(sum(P0_2{j},2)));
% Look for poinst in the overlap
    if (m(1) == 1) 
        list_ic = [];
    else
   if (i==1), 
      list_ic=find(sum(P_1{i},2)+sum(P_1{i+1},2)==2);
   elseif (i==m(1))
      list_ic=find(sum(P_1{i-1},2)+sum(P_1{i},2)==2);
   else
      list_ic=find(sum(P_1{i-1},2)+sum(P_1{i},2)+sum(P_1{i+1},2)==2);
   end
   end
%
    if (m(2) == 1)
        list_jc = [];
    else
   if (j==1),
      list_jc=find(sum(P_2{j},2)+sum(P_2{j+1},2)==2);
   elseif (j==m(2))
      list_jc=find(sum(P_2{j-1},2)+sum(P_2{j},2)==2);
   else
      list_jc=find(sum(P_2{j-1},2)+sum(P_2{j},2)+sum(P_2{j+1},2)==2);
   end
   end
% Build the coarse grid for sub-domain (i,j)
% list_coarse is the list of all coarse points
    for ic=list_ic'
        if (ic == xmin-1) || (ic == xmax+1)
            continue;
        end
       for jc=ymin:stride:ymax
            if (jc == ymin-1) || (jc == ymax+1)
                continue;
            end
         if (coarse_handeled(ic,jc)~=1) %not analyzed 
            coarse_handeled(ic,jc)=1;
            MM = zeros(n,n); 
            [xx,yy]=ndgrid([xmin-1,ic,xmax+1],[ymin-1,jc,ymax+1]);zz=zeros(3,3);zz(2,2)=1;
            F=griddedInterpolant(xx,yy,zz);
            [xx,yy]=ndgrid((xmin-1):(xmax+1),(ymin-1):(ymax+1));zz=F(xx,yy);%mesh(xx,yy,zz);
            zz=zz(2:end-1,2:end-1);MM(xmin:xmax,ymin:ymax)=zz;
            Pc2=[Pc2,MM(:)];
         end
       end
    end
    for jc=list_jc'
        if (jc == ymin-1) || (jc == ymax+1)
            continue;
        end
       for ic=xmin:stride:xmax
           if (ic == xmin-1) || (ic == xmax+1)
               continue;
           end
         if (coarse_handeled(ic,jc)~=1) %not analyzed 
            coarse_handeled(ic,jc)=1;
            MM = zeros(n,n); 
            [xx,yy]=ndgrid([xmin-1,ic,xmax+1],[ymin-1,jc,ymax+1]);zz=zeros(3,3);zz(2,2)=1;
            F=griddedInterpolant(xx,yy,zz);
            [xx,yy]=ndgrid((xmin-1):(xmax+1),(ymin-1):(ymax+1));zz=F(xx,yy);%mesh(xx,yy,zz);
            zz=zz(2:end-1,2:end-1);MM(xmin:xmax,ymin:ymax)=zz;
            Pc2=[Pc2,MM(:)];
         end
       end
    end
%
end
    ind2=ind2+1;
  end
end

Pc2 = sparse(Pc2);

if (coarse == 6 | coarse == 8),
% Pc2=Pc2*randn(size(Pc2,2),min(2*m^2,size(Pc2,2))); not a good idea to introduce randomness in geometry
size(Pc2,2) % number of vectors in the fixed coarse space
%ratio_ff=(size(Pc2,2)/m^2);% number of vectors in the subspace step based subpsace
ratio_ff=(size(Pc2,2)/(m(1)*m(2)));% number of vectors in the subspace step based subpsace
end

