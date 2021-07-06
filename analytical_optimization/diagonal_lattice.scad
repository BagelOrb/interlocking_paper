
L = 3.6;
wb = 1;//2.69;


tpla = true; pp = false; material = tpla;

n=2;



brim = 7;
outer_brim = 20;
side_brim = 30;

tot_l = 50;
rep_z = 5;
rep_y = 5;

specimen_numbers = true;


wmin = .3;
hmin = .2;


wa = 2 * wmin;

w = wa + wb;
d = 2*L*w / sqrt(4*L^2-w^2);
da = wa * d/w;
db = d - da;

tot_w = rep_y*d+da;
tot_h = rep_z*2*hmin;

vec = [L,.5*d];
dvec = vec / sqrt(L^2 + .25*d^2) * .2;


module leg()
{
    linear_extrude(hmin+.0001)
    polygon([
    [0,0],
    [L,.5*d],
    [L,.5*d+da-.2],
    [L,.5*d+da]-dvec,
    [0,da]+dvec,
    [0,da+.2]
    ]);
}

module cell()
{
    leg();
    translate([0,d+da,hmin])
    scale([1,-1,1])
    leg();
}
module interface_a()
{
    for (z = [0:rep_z-1])
        for (y = [0:rep_y-1])
        {
            translate([0,d*y,z*2*hmin]) cell();
        }
    for ( i = [0,1])
        translate([0,i*(tot_w - 2*wmin),0])
        cube([.5*L,2*wmin,rep_z*2*hmin]);
}

module mat_a()
{
    interface_a();
    translate([-tot_l,0,0])
    cube([tot_l+.0001,tot_w,tot_h]);
}

module mat_b()
{
    difference()
    {
        translate([.0001,0,0.0001])
        cube([L+tot_l,tot_w,tot_h-.0001]);
        interface_a();
    }
}

module mat(material)
{
    if (material == tpla) mat_a();
    else mat_b();
}

mat(material);

