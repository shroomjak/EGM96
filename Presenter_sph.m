function Presenter_sph
  input = fopen('potential.txt', 'rt');
  N = fscanf(input, '%d', 1);
  u = zeros(N, N);
  for i = 1:N
    for j = 1:N
      LON(i, j) = fscanf(input, '%e', 1);
      LAT(i, j) = fscanf(input, '%e', 1);
      u(i, j) = fscanf(input, '%e', 1);
    endfor
  endfor
  LAT_raw = LAT(:, 1);
  LON_raw = LON(1, :);
  [LON, LAT] = meshgrid(LON_raw, LAT_raw(end:-1:1));
  R_EQ = 6370000;
  [X, Y, Z] = sph2cart(LON, LAT, R_EQ+u);
  surf(X, Y, Z, u, 'EdgeColor', 'none');
  axis equal;
  colorbar;
  shading interp;
  
  colormap(jet)
  colorbar

  title('Geoid Undulation [m]');
endfunction