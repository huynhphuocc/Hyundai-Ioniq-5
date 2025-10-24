% Sample data
[x, y] = meshgrid(-2:0.05:2);
z = exp(-x.^2 - y.^2).*sin(2*x).*cos(2*y);

figure;
surf(x, y, z, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(parula);                    % hoặc 'turbo' nếu muốn kiểu LS-DYNA
shading interp;
axis equal;
view(45, 25);                        % góc nhìn 3D
set(gcf, 'Color', [0.05 0.05 0.05]); % figure dark gray
ax = gca;
ax.Color = [0.1 0.1 0.1];            % plot region background
ax.XColor = [1 1 1];
ax.YColor = [1 1 1];
ax.ZColor = [1 1 1];
ax.GridColor = [0.5 0.5 0.5];
ax.MinorGridColor = [0.3 0.3 0.3];
ax.BoxStyle = 'full';
xlabel('X [m]', 'Color','w');
ylabel('Y [m]', 'Color','w');
zlabel('Deformation [m]', 'Color','w');
title('3D Surface Visualization (Dark Style)', 'Color','w');

% 💡 Thêm ánh sáng phản chiếu (giống ANSYS shading
camlight('headlight');
lighting phong;       % hoặc 'gouraud' nếu muốn mịn hơn
material dull;        % 'shiny' nếu muốn bóng loáng
                                                   