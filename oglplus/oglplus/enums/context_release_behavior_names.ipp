//  File implement/oglplus/enums/context_release_behavior_names.ipp
//
//  Automatically generated file, DO NOT modify manually.
//  Edit the source 'source/enums/oglplus/context_release_behavior.txt'
//  or the 'source/enums/make_enum.py' script instead.
//
//  Copyright 2010-2017 Matus Chochlik.
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt
//
namespace enums {
OGLPLUS_LIB_FUNC StrCRef ValueName_(
	ContextReleaseBehavior*,
	GLenum value
)
#if (!OGLPLUS_LINK_LIBRARY || defined(OGLPLUS_IMPLEMENTING_LIBRARY)) && \
	!defined(OGLPLUS_IMPL_EVN_CONTEXTRELEASEBEHAVIOR)
#define OGLPLUS_IMPL_EVN_CONTEXTRELEASEBEHAVIOR
{
switch(value)
{
#if defined GL_NONE
	case GL_NONE: return StrCRef("NONE");
#endif
#if defined GL_CONTEXT_RELEASE_BEHAVIOR_FLUSH
	case GL_CONTEXT_RELEASE_BEHAVIOR_FLUSH: return StrCRef("CONTEXT_RELEASE_BEHAVIOR_FLUSH");
#endif
	default:;
}
OGLPLUS_FAKE_USE(value);
return StrCRef();
}
#else
;
#endif
} // namespace enums

