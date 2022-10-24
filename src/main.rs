
extern crate glfw;
extern crate gl;
extern crate glm;

use std::time::Instant;
use std::{thread};
use std::sync::mpsc::{channel, Sender};

pub mod renderer;

use glm::vec3;
use lazy_static::lazy_static;
use renderer::{Renderer, Color, color, WINDOW_WIDTH, WINDOW_HEIGHT};


const SAMPLES_PER_PIXEL : u32 = 5000;       // The total number of times to process each pixel
const SAMPLES_PER_PASS : u32 = 500;          // How many times to process a pixel before presenting to the screen

const NUM_THREADS : usize = 16;              // The number of threads to run on


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

#[derive(Clone, Copy)]
struct PixelRow {
    y: u32,
    c: [Color; WINDOW_WIDTH]
}


lazy_static! {
    static ref SCENE : Scene = Scene{
        spheres: vec![
            sphere(vec3(0.0, -500.0, 0.0), 500.0, vec3(0.6, 0.6, 0.6), 0.0),
            sphere(vec3(0.0, 0.5, -2.0), 0.5, vec3(1.0,0.0,0.0), 0.0), 
            sphere(vec3(-0.6, 0.2, -1.7), 0.2, vec3(1.0,1.0,1.0), 0.0), 
            sphere(vec3(4.0, 4.0, -1.0), 1.0, vec3(1.0, 1.0,1.0), 10.0),
            //sphere(vec3(-3.0, 2.0, -2.0), 1.0, vec3(0.2, 0.2, 1.0), 10.0),

            sphere(vec3(0.0, 1.0, -4.0), 1.0, vec3(0.0, 0.0, 1.0), 0.0),
            sphere(vec3(1.0, 0.25, -1.0), 0.2, vec3(1.0, 1.0, 0.0), 0.0)]
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


impl Camera {
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

                const NUM_PASSES : u32 =  SAMPLES_PER_PIXEL / SAMPLES_PER_PASS;

                for i in 0..NUM_PASSES {
                    for y in start_y..end_y {

                        let mut row = [color(0, 0, 0); WINDOW_WIDTH];

                        for x in 0..WINDOW_WIDTH {
                            let local_pixel_index = (y - start_y) * WINDOW_WIDTH + x;

                            let mut color = local_color_cache[local_pixel_index];

                            for _ in 0..SAMPLES_PER_PASS {
                                let x_offset = fastrand::f32();
                                let y_offset = fastrand::f32();

                                let u = (x as f32 + x_offset) / ((WINDOW_WIDTH - 1) as f32);
                                let v = (y as f32 + y_offset) / ((WINDOW_HEIGHT - 1) as f32);
                    
                                let ray = Ray{origin: pos, dir: lower_left_corner + horizontal*u + vertical*v - pos};
                                
                                color = color + Camera::cast_ray(&ray, 8);
                            }
                            
                            local_color_cache[local_pixel_index] = color;

                            let scale = 1.0 / ((i + 1) * SAMPLES_PER_PASS) as f32;
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


                println!("Thread finished in {:.2?}!", start_time.elapsed());
            });
        }
     
    }
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


    while !renderer.should_close() {
        for _ in 0..NUM_THREADS * 4 {
            if let Ok(j) = rx.recv() {
                for x in 0..WINDOW_WIDTH {
                    renderer.set_pixel(x as u32, j.y, &j.c[x]);
                }
            }
        }

        renderer.update();
    }
}
