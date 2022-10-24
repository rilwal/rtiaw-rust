
extern crate glfw;
extern crate gl;
extern crate glm;

use std::time::Instant;
use std::{thread};
use std::sync::mpsc::{channel, Sender};
use std::sync::atomic::Ordering::Acquire;


pub mod renderer;

use glm::{vec3};
use lazy_static::lazy_static;
use renderer::{Renderer, Color, color, WINDOW_WIDTH, WINDOW_HEIGHT};

const RESTIR_ENABLE : bool = false;
const SAMPLES_PER_PIXEL_RESTIR : u32 = 20;       // The total number of times to process each pixel
const SAMPLES_PER_PASS_RESTIR : u32 = 10;          // How many times to process a pixel before presenting to the screen
const RIS_NUM_PICK : usize = 8;             // The number of samples to use from each reservoir
const RIS_NUM_SAMPLES : usize = 32;         // The number of samples to insert into each reservoir

const SAMPLES_PER_PIXEL : u32 = 5000;       // The total number of times to process each pixel
const SAMPLES_PER_PASS : u32 = 500;          // How many times to process a pixel before presenting to the screen

const NUM_THREADS : usize = 16;              // The number of threads to run on
const FIB_SPHERE_POINTS : usize = 1000;     // The number of points to generate as candidates on each light source


static mut ActiveThreads : std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(NUM_THREADS);

#[derive(Copy, Clone)]
struct Material {
    color: glm::Vec3,
    brightness: f32 // for emission

}

struct Sphere {
    center: glm::Vec3,
    radius: f32,
    mat: Material
}

fn sphere(center: glm::Vec3, radius: f32, color: glm::Vec3, brightness: f32) -> Sphere {
    Sphere {center: center, radius: radius, mat: Material{color: color, brightness: brightness}}
}


struct Scene {
    spheres: Vec<Sphere>
}


struct Ray {
    origin: glm::Vec3,
    dir: glm::Vec3
}


struct Camera {
    pos: glm::Vec3,
    focal_length: f32
}


#[derive(Copy, Clone)]
struct HitRecord {
    dist: f32,
    norm: glm::Vec3,
    point: glm::Vec3,
    mat: Material
}

// Struct for reservoir sampling!
// Consider larger resorvoirs for more performance!
// Investigate the performance impact of Option here!
struct Reservoir<T: Copy> {
    current : Option<T>,
    total_weight : f32,
}

impl<T: Copy> Reservoir<T> {
    fn new() -> Self {
        Self{current: None, total_weight: 0.0}
    }

    fn get(&self) -> &Option<T> {
        &self.current
    }

    fn add(&mut self, item: T, weight: f32) {
        if self.current.is_none() {
            self.current = Some(item);
        }

        //else if fastrand::f32() * (self.total_weight + weight) > self.total_weight {
        else if fastrand::f32() <= weight / self.total_weight {

            self.current = Some(item);
        }

        self.total_weight += weight;
    }
}


struct SizedReservoir<T: Copy, const SIZE : usize> {
    current : [Option<T>; SIZE],
    total_weight : f32,
    count: usize
}

impl<T: Copy, const SIZE : usize> SizedReservoir<T, SIZE> {
    fn new() -> Self {
        Self{current: [None; SIZE], total_weight: 0.0, count: 0}
    }

    fn get(&self) -> &[Option<T>; SIZE] {
        &self.current
    }

    fn add(&mut self, item: T, weight: f32) {
        if self.count < SIZE {
            self.current[self.count] = Some(item);
            self.count += 1;
        }

        else {

            if fastrand::f32() <= weight / self.total_weight {
                self.current[fastrand::usize(0..SIZE)] = Some(item);
            }
        }

        self.total_weight += weight;
    }
}


lazy_static! {
    static ref SCENE : Scene = Scene{
        spheres: vec![
            sphere(vec3(0.0, -500.0, 0.0), 500.0, vec3(0.6, 0.6, 0.6), 0.0),
            sphere(vec3(0.0, 0.5, -2.0), 0.5, vec3(1.0,0.0,0.0), 0.0), 
            sphere(vec3(-0.6, 0.2, -1.7), 0.2, vec3(1.0,1.0,1.0), 0.0), 
            sphere(vec3(4.0, 4.0, -1.0), 1.0, vec3(1.0, 1.0,1.0), 10.0),
            //sphere(vec3(-3.0, 2.0, -2.0), 1.0, vec3(0.2, 0.2, 1.0), 10.0),


            //sphere(vec3(0.0, 10.0, -2.0), 0.001, vec3(1.0, 1.0, 1.0), 100.0),
           // sphere(vec3(0.0, 10.0, -1.0), 1.0, vec3(1.0, 1.0, 1.0), 40.0),
            sphere(vec3(0.0, 1.0, -4.0), 1.0, vec3(0.0, 0.0, 1.0), 0.0),
            sphere(vec3(1.0, 0.25, -1.0), 0.2, vec3(1.0, 1.0, 0.0), 0.0)]
    };

    static ref LIGHTS : Vec<&'static Sphere> = {
        SCENE.spheres.iter().filter(|x| {
            x.mat.brightness > 0.01
        }).collect::<Vec<&Sphere>>()
    };

    static ref FIB_SPHERE : Vec<glm::Vec3> = {
        let mut points = Vec::<glm::Vec3>::with_capacity(FIB_SPHERE_POINTS);
        let phi = std::f32::consts::PI * (3.0 - 5.0_f32.sqrt());
        
        for i in 0..FIB_SPHERE_POINTS {
            let y = 1.0 - (i as f32 / (FIB_SPHERE_POINTS - 1) as f32) * 2.0;
            let radius = (1.0 - y * y).sqrt();
            let theta = phi * i as f32;
    
            let x = theta.cos() * radius;
            let z = theta.sin() * radius;
    
            points.push(vec3(x, y, z));
        }
    
        points   
    };
}


fn random_vec3() -> glm::Vec3 {
    vec3(fastrand::f32() * 2.0 - 1.0, fastrand::f32() * 2.0 - 1.0, fastrand::f32() * 2.0 - 1.0)
}


fn random_vec3_in_unit() -> glm::Vec3 {
    let mut r = random_vec3();

    loop {
        if glm::dot(r, r) >= 1.0 {
            break
        }
        r = random_vec3();
    }

    r
}


fn random_vec3_in_hemi(normal : glm::Vec3) -> glm::Vec3 {
    let r = glm::normalize(random_vec3_in_unit());

    match glm::dot(r, normal) > 0.0 {
        true => r,
        false => -r
    }
}



fn ray_sphere_intersection(ray: &Ray, sphere: &Sphere, t_min : f32, t_max : f32) -> Option<HitRecord> {
    let oc = ray.origin - sphere.center;
    let a = glm::dot(ray.dir, ray.dir);
    let half_b = glm::dot(oc, ray.dir);
    let c = glm::dot(oc, oc) - sphere.radius * sphere.radius;

    let disc = half_b*half_b - a*c;

    if disc < 0.0 {
        return None;
    }

    let sqrtd = disc.sqrt();

    let mut root = (-half_b - sqrtd) / a;
    if root < t_min || t_max < root {
        root = (-half_b + sqrtd) / a;
        if root < t_min || t_max < root {
            return None;
        }
    }


    let dist = root;
    let point = ray.origin + ray.dir * dist;
    let norm = (point - sphere.center) / sphere.radius;
    return Some(HitRecord{dist: dist, norm: norm, point: point, mat: sphere.mat});
    
}


fn random_point_on_light() -> (&'static Sphere, glm::Vec3) {
    // First, select a random light and a random point on a unit sphere
    let light = LIGHTS[fastrand::usize(0..LIGHTS.len())];
    let point = FIB_SPHERE[fastrand::usize(0..FIB_SPHERE.len())];

    //let's try res sampling here
    // let mut res = Reservoir::new();
    // for light in LIGHTS.iter() {
    //     res.add(light, light.radius.powi(2));
    // }
    // let light = res.get().unwrap();

    //let light = LIGHTS[random() % LIGHTS.len()];
    //let point = FIB_SPHERE[random() % FIB_SPHERE.len()];

    // Then, return it as a point on the light
    (light, light.center + point * light.radius)
}


impl Camera {

    // Returns only the distance to hit for a given ray
    pub fn ray_hits(ray: &Ray, target: &Sphere) -> bool {
        let mut closest_hit: Option<HitRecord> = None;
        let mut closest_sphere: Option<&Sphere> = None;

        for sphere in &SCENE.spheres {   
            if let Some(hit) = ray_sphere_intersection(ray, &sphere, 0.001, 10000.0) {
                if closest_hit.is_none() {
                    closest_hit = Some(hit);
                    closest_sphere = Some(sphere);
                } 
                if hit.dist < closest_hit.as_ref().unwrap().dist{
                    closest_hit = Some(hit);
                    closest_sphere = Some(sphere);                
                }
            }
        }

        match closest_sphere {
            None => false,
            Some(x) => std::ptr::eq(target, x)
        }
    }

    pub fn cast_ray(ray: &Ray, bounces: usize) -> glm::Vec3 {
        let mut color = vec3(0.0, 0.0, 0.0);

        if bounces <= 0 {
            return color;
        }

        let mut closest_hit: Option<HitRecord> = None;

        for sphere in &SCENE.spheres {    
            if let Some(hit) = ray_sphere_intersection(ray, &sphere, 0.001, 10000.0) {
                if closest_hit.is_none() {
                    closest_hit = Some(hit)
                } 
                if hit.dist < closest_hit.as_ref().unwrap().dist{
                    closest_hit = Some(hit);
                }
            }
        }

        if let Some(hit) = closest_hit {
            // do some naive stuff here, none of that reSTIR for me!
            let diffuse_reflection = Ray{origin: hit.point, dir: random_vec3_in_hemi(hit.norm)};
            let l_dot_n = glm::dot(hit.norm, diffuse_reflection.dir);

            color = color + hit.mat.color * hit.mat.brightness;
            color = color + hit.mat.color * Camera::cast_ray(&diffuse_reflection, bounces - 1) * l_dot_n;
        }
        
        color
    }


    pub fn cast_ray_shitty(ray: &Ray, bounces: usize) -> glm::Vec3 {
        let mut color = vec3(0.0, 0.0, 0.0);

        if bounces <= 0 {
            return color;
        }

        let mut closest_hit: Option<HitRecord> = None;

        for sphere in &SCENE.spheres {    
            if let Some(hit) = ray_sphere_intersection(ray, &sphere, 0.001, 10000.0) {
                if closest_hit.is_none() {
                    closest_hit = Some(hit)
                } 
                if hit.dist < closest_hit.as_ref().unwrap().dist{
                    closest_hit = Some(hit);
                }
            }
        }

        if let Some(hit) = closest_hit {
            // do some naive stuff here, none of that reSTIR for me!
            //let diffuse_reflection = Ray{origin: hit.point, dir: random_vec3_in_hemi(hit.norm)};
            for light in LIGHTS.iter() {
                let light_vec = glm::normalize(light.center - hit.point);
                
                let shadow_ray = Ray{origin: hit.point, dir: light_vec};
                if Camera::ray_hits(&shadow_ray, light) {
                    let l_dot_n = glm::dot(hit.norm, light_vec);
                    color = color + hit.mat.color * light.mat.color  * l_dot_n;
                }
                
            }
        }
        
        color
    }



    pub fn cast_ray_restir(ray: &Ray) -> glm::Vec3 {
        let mut color = vec3(0.0, 0.0, 0.0);

        let mut closest_hit: Option<HitRecord> = None;

        for sphere in &SCENE.spheres {    
            if let Some(hit) = ray_sphere_intersection(ray, &sphere, 0.001, 10000.0) {
                if closest_hit.is_none() {
                    closest_hit = Some(hit)
                } 
                if hit.dist < closest_hit.as_ref().unwrap().dist{
                    closest_hit = Some(hit);
                }
            }
        }

        if let Some(hit) = closest_hit {
            // Let's try reSTIR type of thing here!

            if hit.mat.brightness > 0.001 {
                return hit.mat.color;
            }

            // For resampled importance sampling, we first take some points on some light sources!
            // let mut samples = Vec::with_capacity(RIS_NUM_SAMPLES);
            
            //let mut res = Reservoir::new();
            let mut res : SizedReservoir<_, RIS_NUM_PICK> = SizedReservoir::new();

            for _ in 0..RIS_NUM_SAMPLES {
                // First, take a random point on a random light
                let (light, point) = random_point_on_light();
                
                // Then calculate our weight
                let light_dir = random_vec3_in_hemi(hit.norm); // glm::normalize(point - hit.point);

                // Cheap weight function
                //let weight = light.mat.brightness * glm::dot(light_dir, hit.norm) * glm::dot(point - hit.point, point - hit.point);

                // Expenny one
                let l_dot_n = glm::dot(light_dir, hit.norm);
                let d_squared = glm::dot(point - hit.point, point - hit.point);
                let unshadowed_light_contrib = hit.mat.color * light.mat.color * light.mat.brightness * (l_dot_n / d_squared);
                let weight = glm::length(unshadowed_light_contrib);
                


               // res.add((c, weight), weight);
                res.add((light, point, unshadowed_light_contrib, weight), weight);

            }            

            let res_output = res.get();

            for i in 0..RIS_NUM_PICK {
                let (light, point, unshadowed_light_contrib, weight) = res_output[i].unwrap();
                

                let light_dir = glm::normalize(point - hit.point);
                let ray = Ray{origin: hit.point, dir: light_dir};

                if Camera::ray_hits(&ray, &light) {
                    //let l_dot_n = glm::dot(light_dir, hit.norm);
                    //let d_squared = glm::dot(point - hit.point, point - hit.point);
                    let w = (res.total_weight / RIS_NUM_SAMPLES as f32) / weight;
                    color = color + unshadowed_light_contrib * w;
                } 
            }

            color = color / RIS_NUM_PICK as f32;
        }
        
        color
    }

    pub fn cast_rays(&self, tx: Sender<PixelRow>) {
        let horizontal = glm::vec3(2.0, 0.0, 0.0);
        let vertical = glm::vec3(0.0, 2.0, 0.0);

        let lower_left_corner = self.pos - horizontal/2.0 - vertical/2.0 - glm::vec3(0.0, 0.0, self.focal_length);

        const ROW_HEIGHT : usize = WINDOW_HEIGHT / NUM_THREADS;


        for thread_id in 0..NUM_THREADS {
            let tx = tx.clone();
            let start_y = ROW_HEIGHT * thread_id;
            let end_y = start_y + ROW_HEIGHT;
            let pos = self.pos;
            
            thread::spawn(move|| {
                let start_time = Instant::now();
                let local_block_size = WINDOW_WIDTH * ROW_HEIGHT;
                let mut local_color_cache: Vec<glm::Vec3> = Vec::with_capacity(local_block_size);
                local_color_cache.resize(local_block_size, glm::vec3(0.0, 0.0, 0.0));

                const NUM_PASSES : u32 = if RESTIR_ENABLE { SAMPLES_PER_PIXEL_RESTIR / SAMPLES_PER_PASS_RESTIR } else { SAMPLES_PER_PIXEL / SAMPLES_PER_PASS };

                for i in 0..NUM_PASSES {
                    for y in start_y..end_y {

                        let mut row = [color(0, 0, 0); WINDOW_WIDTH];

                        for x in 0..WINDOW_WIDTH {
                            let local_pixel_index = (y - start_y) * WINDOW_WIDTH + x;

                            let mut color = local_color_cache[local_pixel_index];

                            for _ in 0.. if RESTIR_ENABLE { SAMPLES_PER_PASS_RESTIR } else { SAMPLES_PER_PASS } {
                                let x_offset = fastrand::f32();
                                let y_offset = fastrand::f32();

                                let u = (x as f32 + x_offset) / ((WINDOW_WIDTH - 1) as f32);
                                let v = (y as f32 + y_offset) / ((WINDOW_HEIGHT - 1) as f32);
                    
                                let ray = Ray{origin: pos, dir: lower_left_corner + horizontal*u + vertical*v - pos};
                                
                                if (RESTIR_ENABLE) {
                                    color = color + Camera::cast_ray_restir(&ray);
                                } else {
                                    color = color + Camera::cast_ray(&ray, 8);
                                }
                            }
                            
                            local_color_cache[local_pixel_index] = color;

                            let scale = 1.0 / ((i + 1) * if RESTIR_ENABLE { SAMPLES_PER_PASS_RESTIR } else { SAMPLES_PER_PASS }) as f32;
                            color = vec3((scale * color.x).sqrt(), (scale * color.y).sqrt(), (scale * color.z).sqrt());
                            

                            let c = Color{r: (color.x * 255.0) as u8, g: (color.y * 255.0) as u8, b: (color.z * 255.0) as u8};
                            row[x] = c;


                        }

                        let mut result = tx.send(PixelRow{y: y as u32,c: row});
                        while result.is_err() {
                            std::thread::sleep(std::time::Duration::from_millis(5));
                            result = tx.send(PixelRow{y: y as u32,c: row});

                            println!("Error sending pixel! Thread {}. Error: {}", thread_id, result.unwrap_err());
                        }
                    }
                }
                unsafe {
                    ActiveThreads.fetch_sub(1, Acquire);
                }

                println!("Thread finished in {:.2?}!", start_time.elapsed());
            });
        }
     
    }
}

#[derive(Clone, Copy)]
struct PixelRow {
    y: u32,
    c: [Color; WINDOW_WIDTH]
}

fn main() {

    let camera = Camera {
        pos: glm::vec3(0.0, 1.0, 1.0),
        focal_length: 1.0
    };

    let mut renderer = Renderer::create();
    renderer.initialize();


    let (tx, rx) = channel::<PixelRow>();

    camera.cast_rays(tx);

    let mut done : bool = false;

    while !renderer.should_close() {
        for _ in 0..NUM_THREADS * 4 {
            if let Ok(j) = rx.recv() {
                for x in 0..WINDOW_WIDTH {
                    renderer.set_pixel(x as u32, j.y, &j.c[x]);
                }
            }
        }

        unsafe {
            if !done && ActiveThreads.load(Acquire) <= 0 {
                done = true;
                println!("DONE");
            }
        }
        renderer.update();
    }
}
