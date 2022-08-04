function [stdev_orig, stdev_list] = calc_angles

pref = [-35.35
-128.13
-137.28
-153.19
-69.5
115.01
-28.14
-105.89
-116.88
-37.19
-56.72
-139.52
-75.04
-34.82
136.42
-1.37
-145.19
-138.2
-146.6
6.69
113.22
-175.76
106.49
150.66
-67.82
-105.35
106.87
-79.45
-47.46
121.05
-126.86
26.45
-106.21
-145.15
-174.33
166.69
159.27
174.31
-29.91
-125.75
-140.59
-88.36
116.21
116.52
-21.16
60.05
153.84
131.94
99.74
-145.34
168.75
-159.67
-78.6
-78.48
-106.14
74.5
179.5
8.93
40.45
-51.83
178.75
-170.75
-156.72
-52.94
-71.94
93.04
30.03
-87.34
55.06
-126.34
17.44
171.69
33.9
-86.35
66.59
101.59
62.03
-96.91
147.1
135.53
];

off = [20.8488
42.8828
332.3491
96.1188
80.6333
84.0446
27.6443
73.9719
15.7603
19.3293
16.9082
12.3492
17.343
23.4226
18.8087
83.7146
68.6259
51.2073
69.6751
310.9647
245.7287
229.3587
169.9152
181.9492
17.5195
59.6443
95.1422
24.2896
21.9748
81.9684
295.3105
157.0304
216.8441
81.7416
91.7801
29.0687
103.6604
85.8149
27.4058
20.7634
294.6314
17.0762
160.4737
71.5246
81.7972
92.5296
308.5381
309.3371
251.7445
79.6229
304.0059
27.5066
293.3311
79.0191
25.5726
96.0013
97.0981
29.3501
161.7053
151.7682
318.3175
313.2872
23.3296
77.2626
226.5137
35.808
26.8903
154.0596
26.1101
16.4374
156.4618
165.5954
82.9812
300.5131
43.6446
32.4938
27.3835
6.4843
33.8866
142.9327
];

pref_adj = zeros(length(pref), 1);
off_adj = zeros(length(off), 1);

PATHOUT = 'Z:\Users\jerry\SurfAnalysis\';    
outfile = [PATHOUT 'LinSumvsCon_ALL_12.20.04.dat'];

fid = fopen(outfile, 'a');

for i=1:length(pref)
   if pref(i) > 360
      pref(i) = mod(pref(i), 360);
   elseif pref(i) < -360
      pref(i) = mod(pref(i), -360);
   end
   
   if pref(i) < 0
      if pref(i) < -180
         pref_adj(i) = mod(pref(i), 180);
      else
         pref_adj(i) = mod(pref(i), -180);
      end
   elseif pref(i) > 0
      if pref(i) > 180
         pref_adj(i) = mod(pref(i), -180);
      else
         pref_adj(i) = mod(pref(i), 180);
      end
   end
   
   test1 = mod(off(i), 360);
   test2 = mod(off(i), -360);
      
   if abs(pref_adj(i)-test1) < abs(pref_adj(i)-test2)
      off_adj(i) = test1;
   else
      off_adj(i) = test2;
   end
   
   line = sprintf(' %3.2f %3.2f\n', pref_adj(i), off_adj(i));
   fprintf(fid, '%s', [line]);
end

fprintf(fid, '\r\n');
fclose(fid);

diff_rad = (pref_adj - off_adj) * (3.14159/180);
sin_diff = sin(diff_rad);
cos_diff = cos(diff_rad);
sum_sin = sum(sin_diff(:));
sum_cos = sum(cos_diff(:));
n = length(diff_rad);
sum_sin = sum_sin/n;
sum_cos = sum_cos/n;
r = sqrt(sum_sin^2+sum_cos^2);
big_R = n*r;
z = big_R^2/n

[n z]

figure
plot(pref_adj, [off_adj], '.');
diff = abs(pref_adj - off_adj);
stdev_orig = std(diff);

%calculate z-value for Rayleigh's test

%include the permutation technique here:
do_permutation = 0;
if do_permutation == 1
    for i=1:10
       rand_ind = randperm(length(off));
       off_boot = zeros(length(off), 1);
       off_randperm = off(rand_ind);
       for j=1:length(pref)
          test1 = mod(off_randperm(j), 360);
          test2 = mod(off_randperm(j), -360);
         
          if abs(pref_adj(j)-test1) < abs(pref_adj(j)-test2)
             off_boot(j) = test1;
         else
             off_boot(j) = test2;
         end
     end
       diff = abs(pref_adj - off_boot);
       stdev_list(i) = std(diff);
    end

    list = find(stdev_list <= stdev_orig);
    sig = length(list)
    p_val = sig/length(stdev_list)
end