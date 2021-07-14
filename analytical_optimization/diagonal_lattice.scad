
sample_points = [
    [[1.2, 3.6], "v"],
    [[1.8, 3.6], "w"],
    [[2.4, 3.6], "x"],
    [[3.0, 3.6], "y"],
    [[3.6, 3.6], "z"],
    ];


tpla = true; pp = false; material = tpla;

n=1;
N = len(sample_points);

wmin = .3;
h = .2;

r = .15;

tot_l = 50;
rep_z = round(5 / (2*h));
rep_y = 5;

brim = 7;
outer_brim = 20;
side_brim = 30;


specimen_numbers = true;


$fn = 6;

wa = 2 * wmin;



module leg(wa,wb,da,db,d,L,h)
{
    w = wa + wb;
    M = L - 2*r;
    d = da + db;
    da = wa * d/w;
    db = wb * d/w;
    ea = (wa - 2*r) * d/w;
    eb = (wb - 2*r) * d/w;
    f = 2*r * d/w;
    s = r * w/d;
    t = sqrt(r^2 - s^2);
    
    translate([r,.5*f,0])
    difference()
    {
        union()
        {
            linear_extrude(h+.0001)
            polygon([
            [-r,-f],
            [-t,-f+s],
            [M+t,.5*d-s],
            [M+r,.5*d],
            [M+r,.5*d+ea],
            [M-t,.5*d+ea+s],
            [t,ea+f-s],
            [-r,ea+f],
            ]);
            
            translate([M,.5*d,0]) cylinder(r=r, h = h);
            translate([M,.5*d+ea,0]) cylinder(r=r, h = h);
        }
        
        translate([0,-f,-.1]) cylinder(r=r, h = 2*h);
        translate([0,ea+f,-.1]) cylinder(r=r, h = 2*h);
    }
}

module cell(wa,wb,da,db,d,L,h)
{
    leg(wa,wb,da,db,d,L,h);
    translate([0,d+da,h])
    scale([1,-1,1])
    leg(wa,wb,da,db,d,L,h);
}
module interface_a(wa,wb,da,db,d,L,h, tot_w,tot_h)
{
    for (z = [0:rep_z-1])
        for (y = [0:rep_y-1])
        {
            translate([0,d*y,z*2*h]) cell(wa,wb,da,db,d,L,h);
        }
    for ( i = [0,1])
        translate([0,i*(tot_w - 2*wmin),0])
        cube([.5*L,2*wmin,tot_h]);
}

module mat_a(wa,wb,da,db,d,L,h, tot_w,tot_h,tot_l)
{
    intersection()
    {
        cube([L+tot_l,tot_w,tot_h-.0001]);
        interface_a(wa,wb,da,db,d,L,h, tot_w,tot_h);
    }
    translate([-tot_l,0,0])
    cube([tot_l+.0001,tot_w,tot_h]);
}

module mat_b(wa,wb,da,db,d,L,h, tot_w,tot_h,tot_l)
{
    difference()
    {
        translate([.0001,0,0.0001])
        cube([L+tot_l,tot_w,tot_h-.0001]);
        interface_a(wa,wb,da,db,d,L,h, tot_w,tot_h);
    }
}

function get_d(M, w) = 2*M*w / sqrt(4*M^2-w^2);

module sample(tag, wb, L)
{
    w = wa + wb;
    M = L - 2*r;
    d = get_d(M, w);
    da = wa * d/w;
    db = wb * d/w;

    tot_w = rep_y*d+da;
    tot_h = rep_z*2*h;
    echo("total width", tot_w);
    echo("total height", tot_h);

    if (material == tpla)
    {
        difference()
        {
            mat_a(wa,wb,da,db,d,L,h, tot_w,tot_h,tot_l);
            translate([-tot_l + 3, tot_w / 2, tot_h-.2])
            linear_extrude(1.0)
            text(tag,halign="left", valign="center",size = min(5,6*.1*(tot_w)));
        }
    }
    else
    {
        difference()
        {
            mat_b(wa,wb,da,db,d,L,h, tot_w,tot_h,tot_l);
            translate([tot_l + L - 3, tot_w / 2, tot_h-.2])
            rotate([0,0,180])
            linear_extrude(1.0)
            text(tag,halign="left", valign="center",size = min(5,6*.1*(tot_w)));
        }
    }
}



module samples(from, till)
{
    echo((from + N - round((n-1) / 5 * N)) % N);
    wb = sample_points[(from + N - round((n-1) / 5 * N)) % N][0][0];
    L = sample_points[(from + N - round((n-1) / 5 * N)) % N][0][1];
    w = wa + wb;
    M = L - 2*r;
    d = get_d(M,w);
    echo("d=",d);
    da = wa * d/w;
    tot_w = rep_y*d+da;
    echo("wb, L", wb, L);
    sample(str(sample_points[(from + N - round((n-1) / 5 * N)) % N][1], n)
        , wb
        , L);
    if (till > from)
    {
        translate([0, brim + tot_w,0])
        samples(from + 1, till);
    }
}

samples(0, N-1);

if (false)
{
    L = 3;
    h = .2;
    wa = 1;
    wb = 2;
    w = wa + wb;
    M = L - 2*r;
    d = get_d(M, w);
    da = wa * d/w;
    db = d - da;
    leg(wa,wb,da,db,d,L,h);
}