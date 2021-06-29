

alternate = 0;

tpla=true;pp=false; material = tpla;
jigsaw=0;suture=1; type = suture;

approx_h = 5;
approx_w = 17;

b_mult = 3;//4.48;

h_min = .5;

wamin = .6;
wbmin = wamin * b_mult;

// jigsaw settings
t = 35; // angle
m = wamin * .5;
n = wbmin * .5;
sect = 1 / cos(t);
ra = ( m + n*sect - n ) / (2 - sect);
rb = ( n - ra*(cos(t)-1) ) / cos(t);

// suture settings
suture_dw = -0.3;
suture_w2 = wamin - 2 * suture_dw;
suture_w1 = wbmin - 2 * suture_dw;
suture_l = 2.4;

tot_l = 50;

brim = 10;

fn=32;

h = h_min * alternate + (1-alternate) * approx_h;

dwa = ra*cos(t);
dwb = rb*cos(t);
dla = ra*sin(t);
dlb = rb*sin(t);

if (type == jigsaw)
{
    echo("min width a:", 2*(dwa+dwb-rb));
    echo("min width b:", 2*(dwa+dwb-ra));
}
else
{
    echo("min width a:", suture_w2 + 2*suture_dw);
    echo("min width b:", suture_w1 + 2*suture_dw);
}
elem_w = (dwa+dwb)*2 * (1-type) + type * (suture_w1+suture_w2+2*suture_dw);
elem_l = (dla+dlb+ra+rb) * (1-type) + type * suture_l;


rep_x = round(approx_w / elem_w);
rep_z = round(approx_h / h) * alternate + (1-alternate);


tot_w = rep_x * elem_w;
tot_h = rep_z * h;
echo("elem width:", elem_w);
echo("elem length:", elem_l);
echo("total width:", tot_w);
echo("total height:", tot_h);

module half()
{
    if (type == jigsaw) half_jigsaw();
        else half_suture();
}
module half_suture()
{
    //translate()
    linear_extrude(h)
    polygon(points=[
    [.001,0],
    [.001,elem_l], 
    [-.5*suture_w2,elem_l],
    [-.5*suture_w2-suture_dw,0],
    [-.5*suture_w2-suture_dw,0]
    ]);
}
module half_jigsaw()
{
    translate([0,elem_l-ra])
    difference()
    {
        union()
        {
            difference()
            {
                cylinder(r=ra, h=h, $fn=fn);
                translate([0.001,-ra,-h]) cube([rb,ra+rb,h*4]);
            }
            translate([-dwa-dwb,-rb-dla-dlb,0]) cube([dwa+dwb+.001,rb+dlb,h]);
        }
        translate([-dwa-dwb,-dla-dlb,-h]) cylinder(r=rb, h=h*4, $fn=fn);
    }
}

module whole()
{
    translate([.5*elem_w,0])
    union()
    {
        half();
        scale([-1,1,1]) half();
    }
}

module whole_alt()
{
    //translate([0,0])
    union()
    {
        translate([elem_w,0,0]) half();
        scale([-1,1,1]) half();
    }
}

module interface(mat = true)
{
    for (z = [0:rep_z-1])
    {
        translate([0,0,z*h])
        for (x = [0:rep_x-1])
        {
            translate([elem_w*x,0,0]) 
            if ((alternate && ((z % 2) == 0) == mat) || (!alternate && pp)) whole();
                else whole_alt();
                
        }
    }
}
full();
module full()
{
    if (material == tpla)
    {
        translate([-brim,0,0]) cube([brim,elem_l,h_min]);
        translate([tot_w,0,0]) cube([brim,elem_l,h_min]);
        translate([-brim,-tot_l-brim,0]) cube([tot_w+2*brim, tot_l+brim, h_min]);
        
        translate([0,-tot_l,0]) cube([tot_w, tot_l+.001, tot_h]);
        interface(tpla);
    }
    else
    {
        translate([-brim,elem_l,0]) cube([tot_w+2*brim, tot_l+brim, h_min]);
        difference()
        {
            translate([0,.001,.001]) cube([tot_w, tot_l + elem_l, tot_h-.002]);
            interface(tpla);
        }
        //translate([0,elem_l,0]) cube([tot_w, tot_l, tot_h]);
        //translate([0,elem_l,0]) scale([1,-1,1]) interface(b);
    }
}