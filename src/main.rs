
extern crate glfw;
extern crate gl;
extern crate glm;

pub mod renderer;

use glm::vec3;
use renderer::{Renderer, Color, color, WINDOW_WIDTH, WINDOW_HEIGHT};


struct Sphere {
    center: glm::Vec3,
    radius: f32,
}

fn sphere(center: glm::Vec3, radius: f32) -> Sphere {
    Sphere {center: center, radius: radius}
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
    point: glm::Vec3
}


fn ray_sphere_intersection(ray: &Ray, sphere: &Sphere) -> Option<HitRecord> {
    let oc = ray.origin - sphere.center;
    let a = glm::dot(ray.dir, ray.dir);
    let half_b = glm::dot(oc, ray.dir);
    let c = glm::dot(oc, oc) - sphere.radius * sphere.radius;

    let disc = half_b*half_b - a*c;

    if disc < 0.0 {
        return None;
    } else {
        let dist = (-half_b - disc.sqrt()) / a;
        let point = ray.origin + ray.dir * dist;
        let norm = (point - sphere.center) / sphere.radius;
        return Some(HitRecord{dist: dist, norm: norm, point: point});
    }
}


impl Camera {

    pub fn cast_ray(&self, ray: &Ray, scene : &Scene) -> Color {
        let light_pos = glm::vec3(0.0, 3.0, 0.0);

        let mut closestHit: Option<HitRecord> = None;

        for sphere in &scene.spheres {    
            if let Some(hit) = ray_sphere_intersection(ray, &sphere) {
                if closestHit.is_none() {
                    closestHit = Some(hit)
                } 
                if hit.dist < closestHit.as_ref().unwrap().dist{
                    closestHit = Some(hit);
                }
            }
        }

        if closestHit.is_some() {
            unsafe {
                let hit = closestHit.unwrap_unchecked();
                let l = glm::normalize(light_pos - hit.point);
                return color((255.0 * glm::dot(l, hit.norm)) as u8, 0, 0) ;
            }
        }

        let unit_direction = glm::normalize(ray.dir);
        let t = 0.5 * (unit_direction.y + 1.0);
        let c = glm::vec3(1.0, 1.0, 1.0) * (1.0 - t) + glm::vec3(0.5, 0.7, 1.0) * t;


        return color ((c.x * 255.0) as u8, (c.y * 255.0) as u8, (c.z * 255.0) as u8);
    }

    pub fn cast_rays(&self, renderer: &mut Renderer, scene : &Scene) {
        let horizontal = glm::vec3(2.0, 0.0, 0.0);
        let vertical = glm::vec3(0.0, 2.0, 0.0);

        let lower_left_corner = self.pos - horizontal/2.0 - vertical/2.0 - glm::vec3(0.0, 0.0, self.focal_length);

        // TODO: Multithread this loop!
        for y in 0..WINDOW_HEIGHT {
            for x in 0..WINDOW_WIDTH {
                let u = x as f32 / ((WINDOW_WIDTH - 1) as f32);
                let v = y as f32 / ((WINDOW_HEIGHT - 1) as f32);

                let ray = Ray{origin: self.pos, dir: lower_left_corner + horizontal*u + vertical*v - self.pos};
                let color = self.cast_ray(&ray, scene);

                renderer.set_pixel(x as u32, y as u32, &color);
            }
        }
    }
}


fn main() {
    let camera = Camera {
        pos: glm::vec3(0.0, 0.0, 0.0),
        focal_length: 1.0
    };


    let scene = Scene{
        spheres: vec![sphere(vec3(0.0, 0.0, -1.0), 0.5), sphere(vec3(0.3, 0.3, -0.5), 0.2)]
    };

    let mut renderer = Renderer::create();

    renderer.initialize();

    camera.cast_rays(&mut renderer, &scene);


    while !renderer.should_close() {
        renderer.update();
    }
}
