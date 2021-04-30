element_width = 5.31357202331391 ;
PLA_widths = [1.1083263946711084, 0.6] ;
PP_widths = [4.713572023313909, 4.2052456286428015] ;
hf =  0.4 ;
hc =  0.2 ;
PLA_vs = [0.6777685262281444, 0.7999999999999998] ;
PP_vs = [3.4022988505747125, 2.8824638471771666] ;





w = element_width;
wa = PLA_widths;
wb = PP_widths;
N = len(wa);
echo(N);

va = PLA_vs;
vb = [ for (i=[0:N-1]) PP_vs[N-1-i]]; // reversed, cause in this file we always count from the PLA side

finger_count = 5;
repetitions = 5;

corner_rounding = .2;

function accumulate(list, c) = 
    c <= 0?
    0
    : list[c-1] + accumulate(list, c - 1);

function l(n) = accumulate(va, n) + accumulate(vb, n);


coords = [
    [0, 0]
    , [0, wa[0]/2 + corner_rounding]
    , [corner_rounding, wa[0]/2]
    , for (i = [0:N-1]) [vb[i] + l(i), wa[i]/2]
    , [l(N) - corner_rounding, wa[N-1]/2]
    , [l(N), wa[N-1]/2 - corner_rounding]
    , [l(N), 0]
    ];

echo(coords);


module ra()
{
    difference()
    {
        polygon(points = coords);
        translate([corner_rounding, wa[0]/2 + corner_rounding]) circle(r=corner_rounding, $fn=16);
    }
    translate([l(N) - corner_rounding, wa[N-1]/2 - corner_rounding]) circle(r=corner_rounding, $fn=16);
}
module finger() { ra(); mirror([0,1,0]) ra(); }

module pattern() {
    linear_extrude(hf) {
        for (i = [0:finger_count-1])
            translate([0,i * w, 0])
                finger();
        }
    for (i = [0:N-1])
        translate([vb[i] + l(i),0,hf])
            cube([va[i],w * (finger_count-1),hc]); 
}

h = (hc+hf);

translate([-20,0,0]) cube([20,w * (finger_count-1),h*repetitions]);
//difference()
{
    //cube([40,w * (finger_count-1),h*repetitions-.01]);
for (i = [0:repetitions-1])
{
    translate([0,0,h*i]) pattern();
}
}

