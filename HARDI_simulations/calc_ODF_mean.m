
function vt=calc_ODF_mean(sq_ODF_matrix)

[m n]=size(sq_ODF_matrix); % sq_ODF_matrix is the matrix of ODF sqrt (functions on the unit Hilbert sphere)
v1=sq_ODF_matrix;
vt=rand(m,1); % initial ODF
vt=sqrt(vt/sum(vt));

%vti=vt;

diff(1)=1;
err=1e-6;
it=1;

while diff(it)>err
    
    phi=zeros(m,1);
    % compute the log map 
    
    for i=1:n
        phi(:,1)=phi(:,1)+1/n*(v1(:,i)-dot(v1(:,i),vt)*vt)./sqrt(1-dot(v1(:,i),vt)^2)*acos(dot(v1(:,i),vt)); 
    end
    
    phi_norm=norm(phi);
    
    vt=cos(phi_norm)*vt+sin(phi_norm)*phi/phi_norm;
    
    it=it+1;
    
    diff(it)=phi_norm;
    
end

%figure,plot(diff) % to show convergence

return


    
    
    
    
    