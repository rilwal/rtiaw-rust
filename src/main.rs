
extern crate glfw;
extern crate gl;
extern crate glm;

use std::thread;
use std::sync::mpsc::{channel, Sender, Receiver};

pub mod renderer;

use glm::vec3;
use renderer::{Renderer, Color, color, WINDOW_WIDTH, WINDOW_HEIGHT};

const SAMPLES_PER_PIXEL : u32 = 250;
const SAMPLES_PER_PASS : u32 = 10;
const NUM_THREADS : usize = 8;

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

lazy_static::lazy_static! {
    static ref SCENE : Scene = Scene{
        spheres: vec![
            sphere(vec3(0.0, -500.0, 0.0), 500.0, vec3(0.9, 0.9, 0.9), 0.0),
            sphere(vec3(0.0, 0.5, -2.0), 0.5, vec3(1.0,0.0,0.0), 0.6), 
            sphere(vec3(0.0, 10.0, 0.0), 4.0, vec3(1.0, 1.0, 0.8), 0.9),
            sphere(vec3(0.0, 0.2, -1.0), 0.2, vec3(1.0, 1.0, 0.0), 0.0)]
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
    let mut r = glm::normalize(random_vec3_in_unit());

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

    pub fn cast_ray(ray: &Ray, depth : u32) -> glm::Vec3 {
        if depth <= 0 {
            return vec3(0.0, 0.0, 0.0);
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
            let target = random_vec3_in_hemi(hit.norm);
            return hit.mat.color * hit.mat.brightness + hit.mat.color * Camera::cast_ray(&Ray{origin: hit.point, dir: target}, depth - 1) * 0.8;
        }
        

        let unit_direction = glm::normalize(ray.dir);
        let t = 0.5 * (unit_direction.y + 1.0);
        return glm::vec3(1.0, 1.0, 1.0) * (1.0 - t) + glm::vec3(0.5, 0.7, 1.0) * t;

    }

    pub fn cast_rays(&self, tx: Sender<Pixel>) {
        let horizontal = glm::vec3(2.0, 0.0, 0.0);
        let vertical = glm::vec3(0.0, 2.0, 0.0);

        let lower_left_corner = self.pos - horizontal/2.0 - vertical/2.0 - glm::vec3(0.0, 0.0, self.focal_length);

        const ROW_HEIGHT : usize = WINDOW_HEIGHT / NUM_THREADS;

        for i in 0..NUM_THREADS {
            let tx = tx.clone();
            let start_y = ROW_HEIGHT * i;
            let end_y = start_y + ROW_HEIGHT;
            let pos = self.pos;
            
            thread::spawn(move|| {
                let local_block_size = WINDOW_WIDTH * ROW_HEIGHT;
                let mut local_color_cache: Vec<glm::Vec3> = Vec::with_capacity(local_block_size);
                local_color_cache.resize(local_block_size, glm::vec3(0.0, 0.0, 0.0));

                const NUM_PASSES : u32 = SAMPLES_PER_PIXEL / SAMPLES_PER_PASS;

                for i in 0..NUM_PASSES {
                    for y in start_y..end_y {
                        for x in 0..WINDOW_WIDTH {
                            let local_pixel_index = (y - start_y) * WINDOW_WIDTH + x;

                            let mut color = local_color_cache[local_pixel_index];

                            for _ in 0..SAMPLES_PER_PASS {
                                let x_offset = fastrand::f32();
                                let y_offset = fastrand::f32();

                                let u = (x as f32 + x_offset) / ((WINDOW_WIDTH - 1) as f32);
                                let v = (y as f32 + y_offset) / ((WINDOW_HEIGHT - 1) as f32);
                    
                                let ray = Ray{origin: pos, dir: lower_left_corner + horizontal*u + vertical*v - pos};
                                
                                color = color + Camera::cast_ray(&ray, 2);
                            }
                            
                            local_color_cache[local_pixel_index] = color;

                            let scale = 1.0 / (i * SAMPLES_PER_PASS) as f32;
                            color = vec3((scale * color.x).sqrt(), (scale * color.y).sqrt(), (scale * color.z).sqrt());
                            
                            let c = Color{r: (color.x * 255.0) as u8, g: (color.y * 255.0) as u8, b: (color.z * 255.0) as u8};
                            tx.send(Pixel{x: x as u32, y: y as u32, c: c});
                        }
                    }
                }
            });
        }
    }
}

#[derive(Clone, Copy)]
struct Pixel {
    x: u32,
    y: u32,
    c: Color
}

fn main() {
    let camera = Camera {
        pos: glm::vec3(0.0, 1.0, 1.0),
        focal_length: 1.0
    };

    let mut renderer = Renderer::create();

    renderer.initialize();

    let (tx, rx) = channel::<Pixel>();
    camera.cast_rays(tx);

    while !renderer.should_close() {
        for _ in 0..WINDOW_WIDTH * 256 {
            if let Ok(j) = rx.recv() {
                renderer.set_pixel(j.x, j.y, &j.c);
            }
        }

        renderer.update();
    }
}
