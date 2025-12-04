function accel_time = extract_acc(mesh, res)

A = res.A;
zloc = mesh.z_acc;
elemZ = mesh.elemZ;

accel_time = zeros(numel(zloc), size(A,2));

for iLoc = 1:numel(zloc)
    node = find(abs(elemZ - zloc(iLoc)) == min(abs(elemZ - zloc(iLoc))),1);
    dof  = 2*(node-1)+1;
    accel_time(iLoc,:) = A(dof,:);
end
