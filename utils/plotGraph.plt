set xtics ('first' 1, 'second' 2, 'third' 3) rotate by 45 right
plot "tsplib_ans.res" u 0:2:xtic(1) title 'error rate' w lines
