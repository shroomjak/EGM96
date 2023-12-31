function Presenter_flat
  input = fopen('potential.txt', 'rt');
  N = fscanf(input, '%d', 1);
  u = zeros(N, N);
  LON = zeros(N, N);
  LAT = zeros(N, N);
  for i = 1:N
    for j = 1:N
      LON(i, j) = fscanf(input, '%e', 1);
      LAT(i, j) = fscanf(input, '%e', 1);
      u(i, j) = fscanf(input, '%e', 1);
    endfor
  endfor
  num_max = fscanf(input, '%d', 1);
  maximums = zeros(2, num_max);
  for i = 1:num_max
    maximums(1, i) = fscanf(input, '%e', 1)+1;
    maximums(2, i) = fscanf(input, '%e', 1)+1;
  endfor
  num_min = fscanf(input, '%d', 1);
  minimums = zeros(2, num_min);
  for i = 1:num_min
    minimums(1, i) = fscanf(input, '%e', 1)+1;
    minimums(2, i) = fscanf(input, '%e', 1)+1;
  endfor
  for i = 1:num_min
    min_lon(i) = LON(minimums(1, i), minimums(2, i));
    min_lat(i) = LAT(minimums(1, i), minimums(2, i));
  endfor
  for i = 1:num_max
    max_lon(i) = LON(maximums(1, i), maximums(2, i));
    max_lat(i) = LAT(maximums(1, i), maximums(2, i));
  endfor
  figure(1);
  contourf(LON*180/pi, LAT*180/pi, u, 100, 'linecolor', 'None');
  xlabel('longitude, grad');
  ylabel('latitude, grad');
  xlim([LON(1,1)*180/pi LON(end, end)*180/pi]);
  ylim([LAT(1,1)*180/pi LAT(end, end)*180/pi]);
  title('Geoid indulation, m');
  colorbar;
  shading interp;
  colormap(jet);
  axis equal;
  figure(2);
  hold on;
  plot(min_lon*180/pi, min_lat*180/pi, 'ob');
  plot(max_lon*180/pi, max_lat*180/pi, 'or');
  xlim([LON(1,1)*180/pi LON(end, end)*180/pi]);
  ylim([LAT(1,1)*180/pi LAT(end, end)*180/pi]);
  xlabel('longitude, grad');
  ylabel('latitude, grad');
  title('Extremums position');
  legend('Minimums', 'Maximums')
  axis equal;
endfunction