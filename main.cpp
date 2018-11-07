#include <glibmm/main.h>
#include <gtkmm/glarea.h>
#include <gtkmm/window.h>


//#include <gtkmm/builder.h>

#include "calc.hpp"
#include "shader.hpp"
#include <vector>


#include <optional>

#include <array>
#include <chrono>
#include <random>
//#include <oglplus/math/matrix.hpp>
using namespace std::complex_literals;


constexpr int field_division_half = 25;
constexpr double field_range_half = 1.5;
constexpr double tick = field_range_half / field_division_half;


int division = 500;


int main()
{
    auto app = Gtk::Application::create();

    Gtk::Window window;
    Gtk::GLArea glarea;
    window.set_default_size(800, 800);

    window.add(glarea);

    std::optional<Shader> shader;
    std::optional<Shader::VertexArray> va;
    std::optional<Shader::VertexArray> va_dots;

    std::vector<Shader::Vertex> buffer;
    std::vector<Shader::Vertex> buffer_dots;


    auto push_dot = [&](const Shader::Vertex& vertex) {
        buffer_dots.push_back(vertex);
    };

    auto push_line = [&](const Shader::Vertex& v1, const Shader::Vertex& v2) {
        buffer.push_back(v1);
        buffer.push_back(v2);
    };


    std::vector<Dipole> dipoles{{0.00 /*1*/, 1}, {0.00 /*1*/, -1}};
    std::vector<Dipole> estimated;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::normal_distribution distrib;

    glarea.signal_render().connect([&](const Glib::RefPtr<Gdk::GLContext>& /*context*/) {
        glarea.make_current();
        glarea.throw_if_error();
        auto viewport = oglplus::Context::Viewport();
        double aspect_ratio = double(viewport.Width()) / viewport.Height();

        oglplus::Context::Clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        oglplus::Context::ClearColor(1.0, 1.0, 1.0, 1.0);
        oglplus::Context::CullFace(oglplus::Face::FrontAndBack);

        oglplus::Mat4f m;
        m.Set(0, 0, float(1 / (field_range_half * aspect_ratio)));
        m.Set(1, 1, float(1 / field_range_half));
        shader->setMatrix(m);
        va->draw();
        va_dots->draw();
        return true;
    });

    glarea.signal_unrealize().connect([&] {
        glarea.make_current();
        glarea.throw_if_error();

        //release opengl objects
        shader = std::nullopt;
        va = std::nullopt;
        va_dots = std::nullopt;
    },
        false);


    auto vertex = [](std::complex<double> z, double r, double g, double b) {
        return Shader::Vertex{oglplus::Vec2f{float(z.real()), float(z.imag())}, oglplus::Vec3f{float(r), float(g), float(b)}};
    };

    double b_lenscale = tick * 1;
    double m_lenscale = 0.1;


    glarea.signal_realize().connect([&] {
        glarea.make_current();
        glarea.throw_if_error();

        shader = Shader{};

        va = shader->createVertexArray(oglplus::PrimitiveType::Lines, 3);
        va_dots = shader->createVertexArray(oglplus::PrimitiveType::Points, 10);
    });

    double theta_tick = 2 * PI / division;
    auto draw = [&] {
        buffer.clear();
        buffer_dots.clear();

        for (int iy = -field_division_half; iy <= field_division_half; iy++) {
            for (int ix = -field_division_half; ix <= field_division_half; ix++) {
                std::complex<double> z = {ix * tick, iy * tick};
                std::complex<double> b = calcField(dipoles, z);
                b /= std::abs(b);
                push_line(vertex(z, 1.0, 0.0, 0.0),
                    vertex(z + b_lenscale * b, 1.0, 1.0, 0.0));
            }
        }


        for (auto& dipole : dipoles) {
            push_line(vertex(dipole.z, 0, 0, 1), vertex(dipole.z + m_lenscale * dipole.m, 0.8, 0.8, 1));
        }

        for (auto& dipole : estimated) {
            push_line(vertex(dipole.z, 0, 0.5, 0), vertex(dipole.z + m_lenscale * dipole.m, 0.5, 1, 0.5));
        }


        for (int i = 0; i < division; i++) {
            std::complex z = std::exp(std::complex<double>{0, i * theta_tick});
            std::complex z2 = std::exp(std::complex<double>{0, (i + 1) * theta_tick});
            push_line(vertex(z, 1.0, 0.0, 0.0), vertex(z2, 1.0, 0.0, 0.0));
            push_dot(vertex(z, 1.0, 0.0, 0.0));
        }

        va->set(buffer.data(), buffer.size());
        va_dots->set(buffer_dots.data(), buffer_dots.size());
        glarea.queue_draw();
    };

    std::chrono::time_point start_time = std::chrono::steady_clock::now();
    Glib::signal_timeout().connect([&] {
        static bool prev_autoplay = false;
        std::chrono::time_point now = std::chrono::steady_clock::now();
        std::chrono::duration t = now - start_time;
        //        dipoles[0].z = 0.5 * std::exp(double(t.count()) * 1.0e-9i);
        //        dipoles[1].z = 0.7 * std::exp(-1.22 * double(t.count()) * 1.0e-9i);

        auto data = potentialTable(dipoles, division);
        //        for(int i=0;i<data.size();i++){
        //            std::cout<<i<<": "<<data[i]<<"\n";
        //        }
        estimated = estimate(data, 11);

        for (int i = 0; i < data.size(); i++) {
            std::cout << i << ":" << data[i].p << std::endl;
        }

        std::cout << "===\n";
        for (Dipole& d : estimated) {
            std::cout << "m: " << d.m << " z: " << d.z << '\n';
        }
        std::cout << "===\n"
                  << std::endl;

        draw();
        return true;
    },
        10);
    glarea.add_events(Gdk::POINTER_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::SCROLL_MASK | Gdk::SMOOTH_SCROLL_MASK);
    glarea.signal_button_press_event().connect([&](GdkEventButton* e) {
        dipoles.push_back(dipoles[0]);
        dipoles.push_back(dipoles[1]);
        dipoles[0].m = {distrib(rng), distrib(rng)};
        return false;
    });

    glarea.signal_motion_notify_event().connect([&](GdkEventMotion* e) {
        static std::complex<double> last_pointer_position = 0;
        if (glarea.get_height() == 0)
            return false;
        std::complex<double> current_pointer_position
            = {((2 * e->x - glarea.get_width()) / (glarea.get_height())) * field_range_half,
                (2 * e->y - glarea.get_height()) / (-glarea.get_height()) * field_range_half};
        //        Vector2d diff = (current_pointer_position - last_pointer_position) / glarea->get_height();
        /*
        last_pointer_position = current_pointer_position;
         */

        dipoles[0].z = current_pointer_position;
        dipoles[1] = dipoles[0].conjugate();
        return true;
    });
    /*glarea->signal_scroll_event().connect([&](GdkEventScroll* e) {
        constexpr double scroll_scale = 0.05;
        if (e->state == 0) {
            distance *= std::exp(scroll_scale * e->delta_y);
            glarea->queue_draw();
        } else if (e->state == GDK_CONTROL_MASK) {
            distance *= std::exp(scroll_scale * e->delta_y);
            fov /= std::exp(scroll_scale * e->delta_y);
            glarea->queue_draw();
        }
        return true;
    });
    window_state->show();
    window_params->show();
    //    window_result->show();
    */
    glarea.show();
    window.show();
    app->run(window);
    return 0;
}