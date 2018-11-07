#pragma once

#include <cmath>

#include <epoxy/gl.h>

#include <oglplus/config/basic.hpp>
#include <oglplus/shader.hpp>

#include <oglplus/buffer.hpp>
#include <oglplus/context.hpp>
#include <oglplus/math/matrix.hpp>
#include <oglplus/math/vector.hpp>
#include <oglplus/program.hpp>
#include <oglplus/uniform.hpp>
#include <oglplus/vertex_array.hpp>
#include <oglplus/vertex_attrib.hpp>

#include <string>
#include <vector>
//#include <eigen3/unsupported/Eigen/OpenGLSupport>

using namespace std::string_view_literals;


struct Shader {
    oglplus::Program program;
    oglplus::VertexArrayAttrib v_pos;
    oglplus::VertexArrayAttrib v_color;

    oglplus::Uniform<oglplus::Mat4f> v_matrix;

    struct Vertex {
        oglplus::Vec2f pos;
        oglplus::Vec3f color;
    };

    Shader()
        : program{[] {
              try {
                  oglplus::Program ret;
                  oglplus::VertexShader vs(R"(
#version 330

in vec2 pos;
in vec3 color;
uniform mat4 matrix;

out vec4 f_color;
void main() {
  gl_Position = matrix*vec4(pos, 0.0, 1.0);
  f_color=vec4(color, 1.0);
}
)");
                  oglplus::FragmentShader fs(R"(
#version 330
in vec4 f_color;

out vec4 outputColor;
void main() {
  outputColor = f_color;
}
)");
                  vs.Compile();
                  fs.Compile();
                  ret.AttachShader(vs);
                  ret.AttachShader(fs);
                  ret.Link();
                  return ret;
              } catch (oglplus::ProgramBuildError& ex) {
                  std::cerr << ex.Log() << std::endl;
                  throw;
              }
          }()},

          v_pos(program, "pos"),
          v_color(program, "color"),
          v_matrix(program, "matrix")

    {
    }

    void setMatrix(const oglplus::Mat4f& mat)
    {
        program.Use();
        v_matrix.Set(mat);
    }
    struct VertexArray {
        Shader* shader;
        oglplus::VertexArray vao;
        oglplus::Buffer buf;
        oglplus::PrimitiveType primitive_type;
        int len;
        int primitive_size;

        VertexArray(Shader& shader, oglplus::PrimitiveType primitive_type = oglplus::PrimitiveType::LineStrip, int primitive_size = 5) : shader(&shader), primitive_type(primitive_type), primitive_size(primitive_size)
        {
            vao.Bind();
            buf.Bind(oglplus::BufferTarget::Array);
            shader.v_pos.Pointer(3, oglplus::DataType::Float, false, sizeof(Vertex), (void*)offsetof(Vertex, pos));
            shader.v_pos.Enable();

            shader.v_color.Pointer(3, oglplus::DataType::Float, false, sizeof(Vertex), (void*)offsetof(Vertex, color));
            shader.v_color.Enable();
        }

        VertexArray(Shader& shader, Vertex* start, int len,
            oglplus::PrimitiveType primitive_type = oglplus::PrimitiveType::LineStrip,
            oglplus::BufferUsage usage = oglplus::BufferUsage::StaticDraw,
            int primitive_size = 5) : VertexArray(shader, primitive_type, primitive_size)
        {
            set(start, len, usage);
        }
        void set(Vertex* start, int len, oglplus::BufferUsage usage = oglplus::BufferUsage::StaticDraw)
        {
            buf.Bind(oglplus::BufferTarget::Array);
            oglplus::Buffer::Data(oglplus::BufferTarget::Array, sizeof(Vertex) * len, start, usage);
            this->len = len;
        }
        void setPartial(int start_index, Vertex* start, int len)
        {
            buf.Bind(oglplus::BufferTarget::Array);
            oglplus::Buffer::SubData(oglplus::BufferTarget::Array, sizeof(Vertex) * start_index, sizeof(Vertex) * len, start);
        }

        void invalidate()
        {
            buf.InvalidateData();
        }
        void invalidatePartial(int start_index, int len)
        {
            buf.InvalidateSubData(sizeof(Vertex) * start_index, sizeof(Vertex) * len);
        }
        void draw()
        {
            draw(0, len);
        }
        void draw(int start_index, int len)
        {
            shader->draw(*this, start_index, len);
        }
    };

    VertexArray createVertexArray(oglplus::PrimitiveType primitive_type = oglplus::PrimitiveType::LineStrip, int primitive_size = 5)
    {
        return VertexArray(*this, primitive_type, primitive_size);
    }

    VertexArray createVertexArray(Vertex* start, int len,
        oglplus::PrimitiveType primitive_type = oglplus::PrimitiveType::LineStrip,
        oglplus::BufferUsage usage = oglplus::BufferUsage::StaticDraw,
        int primitive_size = 5)
    {
        return VertexArray(*this, start, len, primitive_type, usage, primitive_size);
    }
    void draw(const VertexArray& vertex_array, int start_index, int len)
    {
        program.Use();
        if (vertex_array.primitive_type == oglplus::PrimitiveType::Points) {
            oglplus::Context::PointSize(vertex_array.primitive_size);
        } else if (vertex_array.primitive_type == oglplus::PrimitiveType::LineLoop
                   || vertex_array.primitive_type == oglplus::PrimitiveType::Lines
                   || vertex_array.primitive_type == oglplus::PrimitiveType::LinesAdjacency
                   || vertex_array.primitive_type == oglplus::PrimitiveType::LineStrip
                   || vertex_array.primitive_type == oglplus::PrimitiveType::LineStripAdjacency) {
            oglplus::Context::LineWidth(vertex_array.primitive_size);
        }
        vertex_array.vao.Bind();
        oglplus::Context::DrawArrays(vertex_array.primitive_type, start_index, len);
    }

    void draw(const VertexArray& vertex_array)
    {
        draw(vertex_array, 0, vertex_array.len);
    }
};
