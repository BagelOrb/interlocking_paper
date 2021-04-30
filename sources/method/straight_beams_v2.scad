pla_w = 0.35 * 2;
pp_w = pla_w * 4.25;// 0.38 * 2;
tot_w = pla_w + pp_w;

h = 0.2;

rep = 3; // repetitions

d_rep = round(sqrt(2)*rep) + 2;

pattern_l = .5 * sqrt(2) * (tot_w * (d_rep - 1));

echo(pattern_l);

// approx dimensions of total bar
hh_=1;
ww_=15;
ll_=50;


finger_l = rep * tot_w;

n = round(ww_/tot_w); // numbre of fingers
ww = n*tot_w;
layer_count = round(hh_ / h / 2);
hh = layer_count * 2 * h;
ll = ll_ - finger_l / 2;


brim_h = 0.19;
brim_w = 7;

module layer()
{
    // fingers
    for (i = [0:n-1])
        translate([i * tot_w,0,0]) cube([pla_w, finger_l, h]);

    // beams
    for (i = [0:rep-1])
        translate([0,i*tot_w + pla_w,h]) cube([ww, pla_w, h]);
}

module brim_pla(l=finger_l)
{
    hull()
    {
        translate([0,-ll,0]) cylinder(r=brim_w, h = brim_h);
        translate([ww,-ll,0]) cylinder(r=brim_w, h = brim_h);
        translate([-brim_w, -10,0]) cube([ww+brim_w*2, 10, brim_h]);
    }
    translate([-brim_w,0,0]) cube([brim_w, l/2,brim_h]);
    translate([ww,0,0]) cube([brim_w, l/2,brim_h]);
}

module brim_pp(l=finger_l)
{
    translate([ww,l,0]) rotate([0,0,180]) brim_pla(l);
}

module pla()
{
    for (i = [0:layer_count-1])
        translate([0,0,i*2*h]) layer();
    translate([0,-ll,0]) cube([ww,ll,hh]);
}

module pp()
{
    difference()
    {
        translate([0,.0001,0]) cube([ww-.0001,ll + finger_l,hh-.0001]);
        pla();
    }
}


module diag_interlock()
{
    intersection()
    {
        translate([0,-sqrt(2)*pla_w,0])
        rotate([0,0,45])
        for (y = [-4:10])
        {
            translate([y*tot_w, -y*tot_w,0]) cube([pla_w, tot_w * d_rep - pp_w,h]);
            translate([y*tot_w, -y*tot_w,0]) cube([tot_w, pla_w,h]);
            translate([y*tot_w, -y*tot_w,h]) cube([tot_w * d_rep - pp_w,pla_w, h]);
            translate([y*tot_w, -y*tot_w,h]) cube([pla_w, tot_w,h]);
        }
        translate([0,-1,0]) cube([ww,20,hh]);
    }
}

module pla_diag()
{
    for (i = [0:layer_count-1])
        translate([0,0,i*2*h]) diag_interlock();

    translate([0,-ll,0]) cube([ww,ll,hh]);
}


module pp_diag()
{
    difference()
    {
        translate([0,.0001,0]) cube([ww-.0001,ll + pattern_l,hh-.0001]);
        pla_diag();
    }
}

color([1,0,0]) pla_diag();
color([1,0,0]) brim_pla(pattern_l);

// ignore WARNING. Just hit F6 to render, cause F5 doesn't work
*color([0,1,0]) pp_diag();
*color([0,1,0]) brim_pp(pattern_l);


//pla();






