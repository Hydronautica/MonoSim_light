function disp_time = extract_disp(mesh, res)

U = res.U;
zloc = mesh.z_disp;
elemZ = mesh.elemZ;

disp_time = zeros(numel(zloc), size(U,2));

for iLoc = 1:numel(zloc)
    node = find(abs(elemZ - zloc(iLoc)) == min(abs(elemZ - zloc(iLoc))),1);
    dof  = 2*(node-1)+1;
    disp_time(iLoc,:) = U(dof,:);
end
