
extern crate glfw;
extern crate gl;
extern crate glm;

use std::{ffi::c_void, mem::size_of};
use glfw::{Action, Context, Key};

#[derive(Copy, Clone)]
struct Color {
    r: u8,
    g: u8,
    b: u8,
}

const vert_shader_src : &str = "
#version 330 core
layout (location = 0) in vec3 aPos;

void main()
{
    gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);
}
";


const frag_shader_src : &str = "
#version 330 core
out vec4 FragColor;

void main()
{
    FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);
} 
";



fn main() {
    let mut glfw = glfw::init(glfw::FAIL_ON_ERRORS).expect("Failed to init GLFW");

    let (mut window, events) = glfw.create_window(300, 300, "Hello this is window", glfw::WindowMode::Windowed)
    .expect("Failed to create GLFW window.");

    window.set_key_polling(true);
    window.make_current();

    gl::load_with(|s| window.get_proc_address(s) as *const _);

    // Image data
    let mut image_data = vec![Color{r: 0, g: 255, b: 0}; 300*300];
    image_data[302] = Color{r: 255, g: 0, b: 255};


    // And model data
    let mut model_data = [glm::vec3(-0.5, -0.5, 0.0), glm::vec3(0.5, -0.5, 0.0), glm::vec3(0.0, 0.5, 0.0)];

    println!("Model data.len(): {}", model_data.len());
    println!("size of vec3: {}", size_of::<glm::Vec3>());

    // And OpenGL objects
    let mut texture : u32 = 0;
    let mut vao : u32 = 0;
    let mut vbo : u32 = 0;
    
    let mut vert_shader_id : u32 = 0;
    let mut frag_shader_id : u32 = 0;

    let mut program_id : u32 = 0;

    unsafe {
        // Setup a VAO
        gl::BindVertexArray(vao);
        let vptr_model = model_data.as_mut_ptr() as *mut c_void;
        let fptr_model = vptr_model as *mut f32;
        println!("First float is {}", *fptr_model);

        // and the vertex buffer object
        gl::GenBuffers(1, &mut vbo);
        gl::BindBuffer(gl::ARRAY_BUFFER, vbo);
        gl::BufferData(gl::ARRAY_BUFFER, (std::mem::size_of::<glm::Vec3>() * model_data.len()) as isize, vptr_model, gl::STATIC_DRAW);

        // and setup the format of the vertex data
        gl::VertexAttribPointer(0, 3, gl::FLOAT, gl::FALSE, 3 * size_of::<f32>() as i32, 0 as *const c_void);
        gl::EnableVertexAttribArray(0);

        // and unset the VAO
        gl::BindVertexArray(0);

        // and lastly, let's setup some shaders!
        vert_shader_id = gl::CreateShader(gl::VERTEX_SHADER);
        frag_shader_id = gl::CreateShader(gl::FRAGMENT_SHADER);

        let vert_shader_src_ptr = vert_shader_src.as_ptr() as *const i8;
        let vert_shader_len = vert_shader_src.len() as i32;

        let frag_shader_src_ptr = frag_shader_src.as_ptr() as *const i8;
        let frag_shader_len = frag_shader_src.len() as i32;

        gl::ShaderSource(vert_shader_id, 1, &vert_shader_src_ptr as *const *const i8, &vert_shader_len); 
        gl::CompileShader(vert_shader_id);

        gl::ShaderSource(frag_shader_id, 1, &frag_shader_src_ptr as *const *const i8, &frag_shader_len); 
        gl::CompileShader(frag_shader_id);

        // and the shader program
        program_id = gl::CreateProgram();
        gl::AttachShader(program_id, vert_shader_id);
        gl::AttachShader(program_id, frag_shader_id);

        gl::LinkProgram(program_id);
    }

    while !window.should_close() {
        // Render the texture here
        unsafe {
            //gl::BufferData(gl::ARRAY_BUFFER, model_data.len(), )

            let (w, h) = window.get_framebuffer_size();
            gl::Viewport(0, 0, w, h);

            gl::UseProgram(program_id);
            gl::Clear(gl::COLOR_BUFFER_BIT | gl::DEPTH_BUFFER_BIT);
            gl::BindVertexArray(vao); // reload all our settings
            gl::DrawArrays(gl::TRIANGLES, 0, 3);

            window.swap_buffers();
        }


        glfw.poll_events();
        for (_, event) in glfw::flush_messages(&events) {
            match event {
                glfw::WindowEvent::Key(Key::Escape, _, Action::Press, _) => {
                    window.set_should_close(true);
                }
                _ => {}
            }
        }
    }
}
