
extern crate glfw;
extern crate gl;
extern crate glm;

pub mod renderer;

use std::{ffi::c_void, mem::size_of};
use glfw::{Action, Context, Key};
use renderer::{Renderer, Color, color, WINDOW_WIDTH, WINDOW_HEIGHT};


struct Sphere {
    center: glm::Vec3,
    radius: f32,
}

struct Ray {
    origin: glm::Vec3,
    dir: glm::Vec3
}


struct Camera {
    pos: glm::Vec3,
    focal_length: f32
}


static mut t_values : Option<Vec<f32>> = None;


fn ray_sphere_intersection(ray: &Ray, sphere: &Sphere) -> bool {
    let oc = ray.origin - sphere.center;
    let a = glm::dot(ray.dir, ray.dir);
    let b = 2.0 * glm::dot(oc, ray.dir);
    let c = glm::dot(oc, oc) - sphere.radius * sphere.radius;

    let disc = b*b - 4.0*a*c;


    return disc > 0.0;
}


impl Camera {

    pub fn cast_ray(&self, ray: &Ray) -> Color {
        let sphere = Sphere{center: glm::vec3(0.0, 0.0, -2.0), radius: 0.5};

        if ray_sphere_intersection(ray, &sphere) {
            return color(255, 0, 0);
        }

        let unit_direction = glm::normalize(ray.dir);
        let t = 0.5 * (unit_direction.y + 1.0);
        let c = glm::vec3(1.0, 1.0, 1.0) * (1.0 - t) + glm::vec3(0.5, 0.7, 1.0) * t;


        return color ((c.x * 255.0) as u8, (c.y * 255.0) as u8, (c.z * 255.0) as u8);
    }

    pub fn cast_rays(&self, renderer: &mut Renderer) {
        let origin = glm::vec3(0.0, 0.0, 0.0);
        let horizontal = glm::vec3(2.0, 0.0, 0.0);
        let vertical = glm::vec3(0.0, 2.0, 0.0);

        let lower_left_corner = origin - horizontal/2.0 - vertical/2.0 - glm::vec3(0.0, 0.0, self.focal_length);

        for y in 0..renderer::WINDOW_HEIGHT {
            for x in 0..renderer::WINDOW_WIDTH {
                let u = x as f32 / ((renderer::WINDOW_WIDTH - 1) as f32);
                let v = y as f32 / ((renderer::WINDOW_HEIGHT - 1) as f32);

                let ray = Ray{origin: origin, dir: lower_left_corner + horizontal*u + vertical*v - origin};
                let color = self.cast_ray(&ray);

                renderer.set_pixel(x as u32, y as u32, &color);
            }
        }
    }
}


fn main() {
    unsafe {
        t_values = Some(Vec::with_capacity((WINDOW_WIDTH * WINDOW_HEIGHT) as usize));
        t_values.as_mut().unwrap().resize(WINDOW_WIDTH * WINDOW_HEIGHT, 0.0);
    }

    let s = Sphere {
        center: glm::vec3(0.0, 0.0, 0.0),
        radius: 1.0
    };

    let camera = Camera {
        pos: glm::vec3(0.0, 0.0, 0.0),
        focal_length: 1.0
    };


    let mut renderer = Renderer::create();

    renderer.initialize();

    camera.cast_rays(&mut renderer);


    while !renderer.should_close() {
        renderer.update();
    }
}
