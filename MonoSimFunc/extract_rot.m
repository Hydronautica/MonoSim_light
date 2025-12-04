function rot_time = extract_rot(mesh, res)

U = res.U;
zloc = mesh.z_rot;
elemZ = mesh.elemZ;

rot_time = zeros(numel(zloc), size(U,2));

for iLoc = 1:numel(zloc)
    node = find(abs(elemZ - zloc(iLoc)) == min(abs(elemZ - zloc(iLoc))),1);
    dof  = 2*(node-1)+2;
    rot_time(iLoc,:) = U(dof,:);
end
